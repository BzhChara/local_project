import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.preprocessing import StandardScaler

# 定义蛋白质特征提取模型
class ProteinFeatureExtractor(nn.Module):
    def __init__(self, input_dim, hidden_dim, num_layers, embedding_dim):
        super(ProteinFeatureExtractor, self).__init__()
        self.embedding = nn.Linear(input_dim, embedding_dim) 
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_dim, 32)

    def forward(self, x):
        x = self.embedding(x)
        x = F.relu(x)
        x = x.unsqueeze(1)
        lstm_out, _ = self.lstm(x)
        protein_feature = lstm_out[:, -1, :]
        protein = self.fc(protein_feature)
        return protein

# 定义寡糖特征提取模型
class GlycanFeatureExtractor(nn.Module):
    def __init__(self, input_dim, embedding_dim):
        super(GlycanFeatureExtractor, self).__init__()
        self.embedding_dim = embedding_dim
        self.embedding = nn.Embedding(input_dim, embedding_dim)
        self.query = nn.Linear(embedding_dim, embedding_dim)
        self.key = nn.Linear(embedding_dim, embedding_dim)
        self.value = nn.Linear(embedding_dim, embedding_dim)
        self.fc = nn.Linear(embedding_dim, 32)

    def forward(self, x):
        embedding_vector = self.embedding(x)
        Q = self.query(embedding_vector)
        K = self.key(embedding_vector)
        V = self.value(embedding_vector)
        K_T = K.transpose(-2, -1)
        attention_weights = F.softmax(torch.matmul(Q, K_T) / (self.embedding_dim ** 0.5), dim=-1)
        attention_output = torch.matmul(attention_weights, V)
        glycan_feature = attention_output.mean(dim=1)  # 取平均值以匹配全连接层的输入形状
        glycan = self.fc(glycan_feature)
        return glycan

# 定义Transformer模型
class TransformerRegressor(nn.Module):
    def __init__(self, input_dim, num_heads, num_encoder_layers, num_classes):
        super(TransformerRegressor, self).__init__()
        self.embedding = nn.Linear(input_dim, num_classes)
        self.encoder_layer = nn.TransformerEncoderLayer(d_model=num_classes, nhead=num_heads, batch_first=True)
        self.transformer_encoder = nn.TransformerEncoder(self.encoder_layer, num_layers=num_encoder_layers)
        self.fc_out = nn.Linear(num_classes, 1)  # 回归任务的输出是1个值
        self.dropout = nn.Dropout(0.5)

    def forward(self, x):
        x = self.embedding(x)
        x = x.unsqueeze(0)  # 添加序列维度
        x = self.transformer_encoder(x)
        x = x.squeeze(0)  # 移除序列维度
        x = self.dropout(x)
        x = self.fc_out(x)
        return x

# 检查CUDA是否可用
if torch.cuda.is_available():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("CUDA is available. You can train your model on GPU.")
else:
    print("CUDA is not available. You will train your model on CPU.")

# 读取数据
combined_data_Siglec = pd.read_csv("combined_encode_Siglec.csv", header=None)
combined_data_lectin = pd.read_csv("combined_encode_lectin.csv", header=None)
combined_data_ConA = pd.read_csv("combined_encode_ConA.csv", header=None)
combined_data_Galectin = pd.read_csv("combined_encode_Galectin.csv", header=None)
combined_data_train = pd.read_csv("combined_encode_train.csv", header=None)

# 初始化特征提取模型
protein_extractor = ProteinFeatureExtractor(input_dim=2710, hidden_dim=64, num_layers=2, embedding_dim=64).to(device)
glycan_extractor = GlycanFeatureExtractor(input_dim=176, embedding_dim=64).to(device)
    
# 超参数设置
input_dim = 64
num_heads = 2
num_encoder_layers = 2
num_classes = 512  # 嵌入维度

model = TransformerRegressor(input_dim, num_heads, num_encoder_layers, num_classes).to(device)

# 定义损失函数和优化器
criterion = nn.MSELoss().to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.0005)

batch_size = 32

# 训练模型
def train_model(model, epochs):
    model.train()
    for epoch in range(epochs):
        for i in range(0, len(combined_data_train), batch_size):  
            batch_data_train = combined_data_train.iloc[i:i+batch_size].values
            protein_data_train = torch.tensor(batch_data_train[:, :2710], dtype=torch.float32).to(device)
            glycan_data_train = torch.tensor(batch_data_train[:, 2711:2886], dtype=torch.long).to(device)
            values_data_train = torch.tensor(batch_data_train[:, 2886], dtype=torch.float32).unsqueeze(1).to(device)
        
            # 提取蛋白质特征
            protein_features = protein_extractor(protein_data_train)
            # 提取寡糖特征
            glycan_features = glycan_extractor(glycan_data_train)
            # 合并特征向量
            combined_features = torch.cat((protein_features, glycan_features, values_data_train), dim=1)
        
            optimizer.zero_grad()
            outputs = model(combined_features[:, :input_dim])
            loss = criterion(outputs, combined_features[:, -1].unsqueeze(-1))
            loss.backward()
            optimizer.step()
        
        print(f'Epoch {epoch+1}/{epochs}, Loss: {loss.detach().cpu().item()}')

train_model(model, 5)

# Siglec
# 评估模型性能
def evaluate_model(model, mse_criterion, mae_criterion):
    model.eval()  # 将模型设置为评估模式
    total_mse_loss = 0  # 初始化总MSE损失
    total_mae_loss = 0  # 初始化总MAE损失
    with torch.no_grad():  # 在评估期间不计算梯度
        for i in range(0, len(combined_data_Siglec), batch_size):  
            batch_data_test = combined_data_Siglec.iloc[i:i+batch_size].values
            protein_data_test = torch.tensor(batch_data_test[:, :2710], dtype=torch.float32).to(device)
            glycan_data_test = torch.tensor(batch_data_test[:, 2711:2886], dtype=torch.long).to(device)
            values_data_test = torch.tensor(batch_data_test[:, 2886], dtype=torch.float32).unsqueeze(1).to(device)
    
            # 提取蛋白质特征
            protein_features = protein_extractor(protein_data_test)
            # 提取寡糖特征
            glycan_features = glycan_extractor(glycan_data_test)
            # 合并特征向量
            combined_features = torch.cat((protein_features, glycan_features, values_data_test), dim=1)

             # 进行预测
            predictions = model(combined_features[:, :input_dim])
            # 计算MSE损失
            labels = combined_features[:, -1].unsqueeze(-1)
            mse_loss = mse_criterion(predictions, labels)
            total_mse_loss += mse_loss.item()  # 累加MSE损失

            # 计算MAE损失
            mae_loss = mae_criterion(predictions, labels)
            total_mae_loss += mae_loss.item()  # 累加MAE损失

    # 计算平均MSE和MAE损失
    avg_mse_loss = total_mse_loss / len(combined_data_Siglec)
    avg_mae_loss = total_mae_loss / len(combined_data_Siglec)
    return avg_mse_loss, avg_mae_loss

# 定义MSE损失函数
mse_criterion = torch.nn.MSELoss()
# 定义MAE损失函数
mae_criterion = torch.nn.L1Loss()

# 调用评估模型性能函数，并计算MSE和MAE
test_mse, test_mae = evaluate_model(model, mse_criterion, mae_criterion)
print("Siglec:")
print(f'Test MSE: {test_mse}')
print(f'Test MAE: {test_mae}')

# lectin
# 评估模型性能
def evaluate_model(model, mse_criterion, mae_criterion):
    model.eval()  # 将模型设置为评估模式
    total_mse_loss = 0  # 初始化总MSE损失
    total_mae_loss = 0  # 初始化总MAE损失
    with torch.no_grad():  # 在评估期间不计算梯度
        for i in range(0, len(combined_data_lectin), batch_size): 
            batch_data_test = combined_data_lectin.iloc[i:i+batch_size].values
            protein_data_test = torch.tensor(batch_data_test[:, :2710], dtype=torch.float32).to(device)
            glycan_data_test = torch.tensor(batch_data_test[:, 2711:2886], dtype=torch.long).to(device)
            values_data_test = torch.tensor(batch_data_test[:, 2886], dtype=torch.float32).unsqueeze(1).to(device)
    
            # 提取蛋白质特征
            protein_features = protein_extractor(protein_data_test)
            # 提取寡糖特征
            glycan_features = glycan_extractor(glycan_data_test)
            # 合并特征向量
            combined_features = torch.cat((protein_features, glycan_features, values_data_test), dim=1)

             # 进行预测
            predictions = model(combined_features[:, :input_dim])
            # 计算MSE损失
            labels = combined_features[:, -1].unsqueeze(-1)
            mse_loss = mse_criterion(predictions, labels)
            total_mse_loss += mse_loss.item()  # 累加MSE损失

            # 计算MAE损失
            mae_loss = mae_criterion(predictions, labels)
            total_mae_loss += mae_loss.item()  # 累加MAE损失

    # 计算平均MSE和MAE损失
    avg_mse_loss = total_mse_loss / len(combined_data_lectin)
    avg_mae_loss = total_mae_loss / len(combined_data_lectin)
    return avg_mse_loss, avg_mae_loss

# 定义MSE损失函数
mse_criterion = torch.nn.MSELoss()
# 定义MAE损失函数
mae_criterion = torch.nn.L1Loss()

# 调用评估模型性能函数，并计算MSE和MAE
test_mse, test_mae = evaluate_model(model, mse_criterion, mae_criterion)
print("C-type lectin:")
print(f'Test MSE: {test_mse}')
print(f'Test MAE: {test_mae}')

# ConA
# 评估模型性能
def evaluate_model(model, mse_criterion, mae_criterion):
    model.eval()  # 将模型设置为评估模式
    total_mse_loss = 0  # 初始化总MSE损失
    total_mae_loss = 0  # 初始化总MAE损失
    with torch.no_grad():  # 在评估期间不计算梯度
        for i in range(0, len(combined_data_ConA), batch_size):  
            batch_data_test = combined_data_ConA.iloc[i:i+batch_size].values
            protein_data_test = torch.tensor(batch_data_test[:, :2710], dtype=torch.float32).to(device)
            glycan_data_test = torch.tensor(batch_data_test[:, 2711:2886], dtype=torch.long).to(device)
            values_data_test = torch.tensor(batch_data_test[:, 2886], dtype=torch.float32).unsqueeze(1).to(device)
    
            # 提取蛋白质特征
            protein_features = protein_extractor(protein_data_test)
            # 提取寡糖特征
            glycan_features = glycan_extractor(glycan_data_test)
            # 合并特征向量
            combined_features = torch.cat((protein_features, glycan_features, values_data_test), dim=1)

             # 进行预测
            predictions = model(combined_features[:, :input_dim])
            # 计算MSE损失
            labels = combined_features[:, -1].unsqueeze(-1)
            mse_loss = mse_criterion(predictions, labels)
            total_mse_loss += mse_loss.item()  # 累加MSE损失

            # 计算MAE损失
            mae_loss = mae_criterion(predictions, labels)
            total_mae_loss += mae_loss.item()  # 累加MAE损失

    # 计算平均MSE和MAE损失
    avg_mse_loss = total_mse_loss / len(combined_data_ConA)
    avg_mae_loss = total_mae_loss / len(combined_data_ConA)
    return avg_mse_loss, avg_mae_loss

# 定义MSE损失函数
mse_criterion = torch.nn.MSELoss()
# 定义MAE损失函数
mae_criterion = torch.nn.L1Loss()

# 调用评估模型性能函数，并计算MSE和MAE
test_mse, test_mae = evaluate_model(model, mse_criterion, mae_criterion)
print("ConA:")
print(f'Test MSE: {test_mse}')
print(f'Test MAE: {test_mae}')

# Galectin
# 评估模型性能
def evaluate_model(model, mse_criterion, mae_criterion):
    model.eval()  # 将模型设置为评估模式
    total_mse_loss = 0  # 初始化总MSE损失
    total_mae_loss = 0  # 初始化总MAE损失
    with torch.no_grad():  # 在评估期间不计算梯度
        for i in range(0, len(combined_data_Galectin), batch_size):  
            batch_data_test = combined_data_Galectin.iloc[i:i+batch_size].values
            protein_data_test = torch.tensor(batch_data_test[:, :2710], dtype=torch.float32).to(device)
            glycan_data_test = torch.tensor(batch_data_test[:, 2711:2886], dtype=torch.long).to(device)
            values_data_test = torch.tensor(batch_data_test[:, 2886], dtype=torch.float32).unsqueeze(1).to(device)
    
            # 提取蛋白质特征
            protein_features = protein_extractor(protein_data_test)
            # 提取寡糖特征
            glycan_features = glycan_extractor(glycan_data_test)
            # 合并特征向量
            combined_features = torch.cat((protein_features, glycan_features, values_data_test), dim=1)

             # 进行预测
            predictions = model(combined_features[:, :input_dim])
            # 计算MSE损失
            labels = combined_features[:, -1].unsqueeze(-1)
            mse_loss = mse_criterion(predictions, labels)
            total_mse_loss += mse_loss.item()  # 累加MSE损失

            # 计算MAE损失
            mae_loss = mae_criterion(predictions, labels)
            total_mae_loss += mae_loss.item()  # 累加MAE损失

    # 计算平均MSE和MAE损失
    avg_mse_loss = total_mse_loss / len(combined_data_Galectin)
    avg_mae_loss = total_mae_loss / len(combined_data_Galectin)
    return avg_mse_loss, avg_mae_loss

# 定义MSE损失函数
mse_criterion = torch.nn.MSELoss()
# 定义MAE损失函数
mae_criterion = torch.nn.L1Loss()

# 调用评估模型性能函数，并计算MSE和MAE
test_mse, test_mae = evaluate_model(model, mse_criterion, mae_criterion)
print("Galectin:")
print(f'Test MSE: {test_mse}')
print(f'Test MAE: {test_mae}')