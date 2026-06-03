import numpy as np
import pandas as pd 
import math
import torch
import torch.nn as nn
import torch.nn.functional as F
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error, r2_score

# 定义蛋白质特征提取模型
class ProteinFeatureExtractor(nn.Module):
    def __init__(self, input_dim, hidden_dim, num_layers, embedding_dim):
        super(ProteinFeatureExtractor, self).__init__()
        self.embedding = nn.Linear(input_dim, embedding_dim) 
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_dim, 32)

    def forward(self, x):
        x = x.unsqueeze(1)
        x = self.embedding(x)
        x = F.relu(x)
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

"""""
# 定义位置编码模块
class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len):
        super(PositionalEncoding, self).__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        self.register_buffer('pe', pe)

    def forward(self, x):
        x = x + self.pe[:, :x.size(1), :]
        return x
"""
        
# 定义Transformer模型
class TransformerRegressor(nn.Module):
    def __init__(self, input_dim, num_heads, num_encoder_layers, num_classes):
        super(TransformerRegressor, self).__init__()
        self.embedding = nn.Linear(input_dim, num_classes)
        # self.positional_encoding = PositionalEncoding(num_classes, max_len=max_len)
        self.encoder_layer = nn.TransformerEncoderLayer(d_model=num_classes, nhead=num_heads, batch_first=True)
        self.transformer_encoder = nn.TransformerEncoder(self.encoder_layer, num_layers=num_encoder_layers)
        self.norm = nn.LayerNorm(num_classes)
        self.dropout = nn.Dropout(0.3)
        self.fc_out = nn.Linear(num_classes, 1)

    def forward(self, x):
        x = self.embedding(x)  # 嵌入层
        # x = self.positional_encoding(x)  # 添加位置编码
        x = self.transformer_encoder(x)  # Transformer 编码器
        x = x[:, -1, :]  # 取最后一个时间步的输出
        x = self.norm(x)
        x = self.dropout(x)  # Dropout
        x = self.fc_out(x)  # 输出层
        return x

# 检查CUDA是否可用
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# 检查Flash Attention是否可用，尝试手动开启
print(torch.cuda.get_device_properties("cuda").name)
torch.backends.cuda.enable_flash_sdp(True)

# 读取数据
combined_data = pd.read_csv("combined_encode.csv", header=None)

# 在读取数据后添加随机打乱
combined_data = combined_data.sample(frac=1).reset_index(drop=True)

# 计算训练集和测试集的大小
train_size = int(0.8 * len(combined_data))
# print(train_size)
test_size = len(combined_data) - train_size
# 打乱数据后固定划分
train_data = combined_data.iloc[:train_size]
test_data = combined_data.iloc[train_size:]

# 初始化特征提取模型
protein_extractor = ProteinFeatureExtractor(input_dim=2710, hidden_dim=128, num_layers=2, embedding_dim=64).to(device)
glycan_extractor = GlycanFeatureExtractor(input_dim=342, embedding_dim=64).to(device)
    
# 初始化模型
input_dim = 32
num_heads = 4
num_encoder_layers = 4
num_classes = 256  
max_len = 2  # 时间步数目
batch_size = 32

model = TransformerRegressor(input_dim, num_heads, num_encoder_layers, num_classes).to(device)

# 定义损失函数和优化器
criterion = nn.MSELoss().to(device)
all_parameters = list(model.parameters()) + list(protein_extractor.parameters()) + list(glycan_extractor.parameters())
optimizer = torch.optim.Adam(all_parameters, lr=0.0001, weight_decay=1e-4)

# 定义学习率调度器
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, factor=0.5, patience=3)

# 训练模型
def train_model(model, epochs):
    model.train()
    train_losses = []
    test_losses = []
    final_true = []
    final_preds = []
    for epoch in range(epochs):
        epoch_loss = 0
        num_train_batches = (train_size + batch_size - 1) // batch_size  # 向上取整
        for i in range(0, len(train_data), batch_size):
            end = min(i + batch_size, train_size)
            batch_data = train_data.iloc[i:end].values
            protein_data = torch.tensor(batch_data[:, :2710], dtype=torch.float32).to(device)
            glycan_data = torch.tensor(batch_data[:, 2710:3052], dtype=torch.long).to(device)
            values_data = torch.tensor(batch_data[:, 3052], dtype=torch.float32).unsqueeze(1).to(device)
        
            protein_features = protein_extractor(protein_data)
            # print("Protein features shape:", protein_features.shape)
            # print("Protein features sample:", protein_features[0, :5])  # 打印前5个特征值
            glycan_features = glycan_extractor(glycan_data)
            # print("Glycan features shape:", glycan_features.shape)
            # print("Glycan features sample:", glycan_features[0, :5])  # 打印前5个特征值
            combined_features = torch.stack((protein_features, glycan_features), dim=1)  
        
            optimizer.zero_grad()
            outputs = model(combined_features)
            loss = criterion(outputs, values_data)
            loss.backward()
            optimizer.step()
            
            epoch_loss += loss.item()
        
        epoch_loss /= num_train_batches  # 按实际批次数量平均
        train_losses.append(epoch_loss)
        
        test_loss, epoch_preds, epoch_true = evaluate_model(model, test_size, batch_size)
        test_losses.append(test_loss)
        
        # 仅在最后一个 epoch 保存结果
        if epoch == epochs - 1:
            final_true = epoch_true
            final_preds = epoch_preds

        scheduler.step(test_loss)
        
        print(f'Epoch {epoch+1}/{epochs}, Train Loss: {epoch_loss}, Test Loss: {test_loss}')
    
    # 绘制最后一次 epoch 的残差图
    plot_residuals(final_true, final_preds)
    plot_learning_curve(train_losses, test_losses)
    
    return train_losses, test_losses

# 评估模型
def evaluate_model(model, test_size, batch_size):
    model.eval()
    test_loss = 0
    predictions = []
    true_values = []
    num_test_batches = (test_size + batch_size - 1) // batch_size  # 向上取整
    with torch.no_grad():
        for i in range(0, len(test_data), batch_size):
            end = min(i + batch_size, len(combined_data))
            batch_data = test_data.iloc[i:end].values
            protein_data = torch.tensor(batch_data[:, :2710], dtype=torch.float32).to(device)
            glycan_data = torch.tensor(batch_data[:, 2710:3052], dtype=torch.long).to(device)
            values_data = torch.tensor(batch_data[:, 3052], dtype=torch.float32).unsqueeze(1).to(device)
        
            protein_features = protein_extractor(protein_data)
            glycan_features = glycan_extractor(glycan_data)
            combined_features = torch.stack((protein_features, glycan_features), dim=1)  
        
            outputs = model(combined_features)
            loss = criterion(outputs, values_data)
            test_loss += loss.item()
            
            predictions.extend(outputs.cpu().numpy().flatten())
            true_values.extend(values_data.cpu().numpy().flatten())
    
    test_loss /= num_test_batches  # 按实际批次数量平均
    
    mse = test_loss
    mae = mean_absolute_error(true_values, predictions)
    r2 = r2_score(true_values, predictions)
    
    print(f'Test MSE: {mse:.4f}, Test MAE: {mae:.4f}, Test R²: {r2:.4f}')
    
    model.train()
    return test_loss, predictions, true_values

# 绘制学习曲线
def plot_learning_curve(train_losses, test_losses):
    plt.figure(figsize=(10, 6))
    plt.plot(train_losses, label='Training Loss')
    plt.plot(test_losses, label='Validation Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Learning Curve')
    plt.legend()
    plt.show()

def plot_residuals(true_values, predictions):
    residuals = np.array(true_values) - np.array(predictions)
    
    plt.figure(figsize=(12, 6))
    
    # 散点图
    plt.subplot(1, 2, 1)
    plt.scatter(predictions, residuals, alpha=0.5, color='blue')
    plt.axhline(y=0, color='red', linestyle='--', linewidth=2)
    plt.xlabel('Predictions')
    plt.ylabel('Residuals')
    plt.title('Residual Plot')
    plt.grid(True)
    
    # 残差分布直方图
    plt.subplot(1, 2, 2)
    plt.hist(residuals, bins=20, color='green', alpha=0.7)
    plt.xlabel('Residuals')
    plt.ylabel('Frequency')
    plt.title('Residual Distribution')
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()

# 训练模型
train_losses, test_losses = train_model(model, 5)