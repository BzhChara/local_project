import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

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
combined_data = pd.read_csv("combined_encode.csv", header=None)

# 计算训练集的大小
train_size = int(0.8 * len(combined_data))

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
        for i in range(0, train_size, batch_size):  # 只遍历前80%的数据
            batch_data = combined_data.iloc[i:i+batch_size].values
            protein_data = torch.tensor(batch_data[:, :2710], dtype=torch.float32).to(device)
            glycan_data = torch.tensor(batch_data[:, 2711:2886], dtype=torch.long).to(device)
            values_data = torch.tensor(batch_data[:, 2886], dtype=torch.float32).unsqueeze(1).to(device)
        
            # 提取蛋白质特征
            protein_features = protein_extractor(protein_data)
            # 提取寡糖特征
            glycan_features = glycan_extractor(glycan_data)
            # 合并特征向量
            combined_features = torch.cat((protein_features, glycan_features, values_data), dim=1)
        
            optimizer.zero_grad()
            outputs = model(combined_features[:, :input_dim])
            loss = criterion(outputs, combined_features[:, -1].unsqueeze(-1))
            loss.backward()
            optimizer.step()
        
        print(f'Epoch {epoch+1}/{epochs}, Loss: {loss.detach().cpu().item()}')

train_model(model, 5)

# 评估模型性能
def evaluate_model(model, mse_criterion, mae_criterion):
    model.eval()  # 将模型设置为评估模式
    total_mse_loss = 0  # 初始化总MSE损失
    total_mae_loss = 0  # 初始化总MAE损失
    with torch.no_grad():  # 在评估期间不计算梯度
        for i in range(train_size, len(combined_data), batch_size):  # 遍历后20%的数据
            batch_data = combined_data.iloc[i:i+batch_size].values
            protein_data = torch.tensor(batch_data[:, :2710], dtype=torch.float32).to(device)
            glycan_data = torch.tensor(batch_data[:, 2711:2886], dtype=torch.long).to(device)
            values_data = torch.tensor(batch_data[:, 2886], dtype=torch.float32).unsqueeze(1).to(device)
    
            # 提取蛋白质特征
            protein_features = protein_extractor(protein_data)
            # 提取寡糖特征
            glycan_features = glycan_extractor(glycan_data)
            # 合并特征向量
            combined_features = torch.cat((protein_features, glycan_features, values_data), dim=1)

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
    avg_mse_loss = total_mse_loss / (len(combined_data) - train_size)
    avg_mae_loss = total_mae_loss / (len(combined_data) - train_size)
    return avg_mse_loss, avg_mae_loss

# 定义MSE损失函数
mse_criterion = torch.nn.MSELoss()
# 定义MAE损失函数
mae_criterion = torch.nn.L1Loss()

# 调用评估模型性能函数，并计算MSE和MAE
test_mse, test_mae = evaluate_model(model, mse_criterion, mae_criterion)
print(f'Test MSE: {test_mse}')
print(f'Test MAE: {test_mae}')