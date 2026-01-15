import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader, random_split
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

# 定义蛋白质特征提取模型
class ProteinFeatureExtractor(nn.Module):
    def __init__(self, input_dim, hidden_dim, num_layers, embedding_dim):
        super(ProteinFeatureExtractor, self).__init__()
        self.embedding = nn.Linear(input_dim, embedding_dim) 
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_dim, 64)

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
        self.fc = nn.Linear(embedding_dim, 64)

    def forward(self, x):
        embedding_vector = self.embedding(x)
        Q = self.query(embedding_vector)
        K = self.key(embedding_vector)
        V = self.value(embedding_vector)
        K_T = K.transpose(-2, -1)
        attention_weights = F.softmax(torch.matmul(Q, K_T) / (self.embedding_dim ** 0.5), dim=-1)
        attention_output = torch.matmul(attention_weights, V)
        glycan_feature = attention_output.mean(dim=1)
        glycan = self.fc(glycan_feature)
        return glycan

# 定义Transformer回归模型
class TransformerRegressionModel(nn.Module):
    def __init__(self, protein_input_dim, glycan_input_dim, hidden_dim=64, num_layers=2, embedding_dim=64):
        super(TransformerRegressionModel, self).__init__()
        self.protein_extractor = ProteinFeatureExtractor(protein_input_dim, hidden_dim, num_layers, embedding_dim)
        self.glycan_extractor = GlycanFeatureExtractor(glycan_input_dim, embedding_dim)
        self.transformer = nn.TransformerEncoder(
            encoder_layer=nn.TransformerEncoderLayer(d_model=64, nhead=8),
            num_layers=2
        )
        self.fc = nn.Linear(64, 1)

    def forward(self, protein_data, glycan_data):
        protein_features = self.protein_extractor(protein_data)
        glycan_features = self.glycan_extractor(glycan_data)
        combined_features = torch.cat((protein_features.unsqueeze(1), glycan_features.unsqueeze(1)), dim=1)
        transformer_output = self.transformer(combined_features)
        output = self.fc(transformer_output.mean(dim=1))
        return output.squeeze()

# 自定义数据集
class CustomDataset(Dataset):
    def __init__(self, data):
        self.data = data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        protein = self.data.iloc[idx, :2710].values.astype(np.float32)
        glycan = self.data.iloc[idx, 2710:2710+342].values.astype(np.int64)
        target = self.data.iloc[idx, -1].astype(np.float32)
        return torch.tensor(protein), torch.tensor(glycan), torch.tensor(target)

# 训练函数
def train_model(model, train_loader, criterion, optimizer, device):
    model.train()
    for batch_idx, (protein, glycan, target) in enumerate(train_loader):
        protein, glycan, target = protein.to(device), glycan.to(device), target.to(device)
        optimizer.zero_grad()
        output = model(protein, glycan)
        loss = criterion(output, target)
        loss.backward()
        optimizer.step()

# 评估函数
def evaluate_model(model, test_loader, criterion, device):
    model.eval()
    predictions = []
    true_values = []
    total_loss = 0
    with torch.no_grad():
        for protein, glycan, target in test_loader:
            protein, glycan, target = protein.to(device), glycan.to(device), target.to(device)
            output = model(protein, glycan)
            loss = criterion(output, target)
            total_loss += loss.item() * protein.size(0)
            predictions.extend(output.cpu().numpy())
            true_values.extend(target.cpu().numpy())
    avg_loss = total_loss / len(test_loader.dataset)
    return avg_loss, np.array(predictions), np.array(true_values)

# 主函数
def main():
    # 读取数据
    combined_data = pd.read_csv("combined_encode.csv", header=None)
    
    # 数据集划分
    dataset = CustomDataset(combined_data)
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    train_dataset, test_dataset = random_split(dataset, [train_size, test_size])
    
    # 创建数据加载器
    batch_size = 32
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    
    # 初始化模型、损失函数和优化器
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = TransformerRegressionModel(protein_input_dim=2710, glycan_input_dim=342).to(device)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    
    # 训练模型
    num_epochs = 50
    for epoch in range(num_epochs):
        train_model(model, train_loader, criterion, optimizer, device)
        train_loss, _, _ = evaluate_model(model, train_loader, criterion, device)
        test_loss, _, _ = evaluate_model(model, test_loader, criterion, device)
        print(f"Epoch {epoch+1}/{num_epochs}, Train Loss: {train_loss:.4f}, Test Loss: {test_loss:.4f}")
    
    # 最终评估
    _, predictions, true_values = evaluate_model(model, test_loader, criterion, device)
    
    # 计算指标
    mae = mean_absolute_error(true_values, predictions)
    mse = mean_squared_error(true_values, predictions)
    r2 = r2_score(true_values, predictions)
    print(f"MAE: {mae:.4f}, MSE: {mse:.4f}, R2: {r2:.4f}")
    
    # 绘制残差图
    residuals = true_values - predictions
    plt.scatter(predictions, residuals, alpha=0.5)
    plt.axhline(y=0, color='r', linestyle='--')
    plt.xlabel('Predicted Values')
    plt.ylabel('Residuals')
    plt.title('Residual Plot')
    plt.show()

if __name__ == "__main__":
    main()