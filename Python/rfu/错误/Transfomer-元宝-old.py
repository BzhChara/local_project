import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader, random_split
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import matplotlib.pyplot as plt

# 设备配置
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# 数据加载与预处理
def load_data():
    combined_data = pd.read_csv("combined_encode.csv", header=None)
    protein_data = combined_data.iloc[:, :2710].values
    glycan_data = combined_data.iloc[:, 2710:2710+342].values
    target = combined_data.iloc[:, -1].values
    
    # 转换为张量
    protein_tensor = torch.LongTensor(protein_data)
    glycan_tensor = torch.FloatTensor(glycan_data)
    target_tensor = torch.FloatTensor(target)
    
    return protein_tensor, glycan_tensor, target_tensor

# 自定义数据集
class CustomDataset(Dataset):
    def __init__(self, protein, glycan, target):
        self.protein = protein
        self.glycan = glycan
        self.target = target

    def __len__(self):
        return len(self.target)

    def __getitem__(self, idx):
        return self.protein[idx], self.glycan[idx], self.target[idx]

# 模型定义
class ProteinFeatureExtractor(nn.Module):
    def __init__(self, vocab_size=21, embedding_dim=64, hidden_dim=128, num_layers=2):
        super().__init__()
        self.embedding = nn.Embedding(vocab_size, embedding_dim)
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_dim, 64)

    def forward(self, x):
        x = self.embedding(x)
        lstm_out, _ = self.lstm(x)
        return self.fc(lstm_out[:, -1, :])

class GlycanFeatureExtractor(nn.Module):
    def __init__(self, input_dim=342, embedding_dim=64):
        super().__init__()
        self.fc = nn.Linear(input_dim, embedding_dim)
        self.query = nn.Linear(embedding_dim, embedding_dim)
        self.key = nn.Linear(embedding_dim, embedding_dim)
        self.value = nn.Linear(embedding_dim, embedding_dim)
        self.fc_out = nn.Linear(embedding_dim, 64)

    def forward(self, x):
        x = self.fc(x)
        x = x.unsqueeze(1)
        Q = self.query(x)
        K = self.key(x)
        V = self.value(x)
        attn = F.softmax(torch.matmul(Q, K.transpose(1,2)) / np.sqrt(64), dim=-1)
        out = torch.matmul(attn, V).squeeze(1)
        return self.fc_out(out)

class TransformerRegression(nn.Module):
    def __init__(self, protein_extractor, glycan_extractor, d_model=128, nhead=8, num_layers=3):
        super().__init__()
        self.protein_extractor = protein_extractor
        self.glycan_extractor = glycan_extractor
        
        self.fusion_fc = nn.Linear(64 * 2, d_model)
        encoder_layer = nn.TransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=512, batch_first=True)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        self.regressor = nn.Sequential(
            nn.Linear(d_model, 64),
            nn.ReLU(),
            nn.Linear(64, 1)
        )

    def forward(self, protein, glycan):
        protein_feat = self.protein_extractor(protein)
        glycan_feat = self.glycan_extractor(glycan)
        combined = torch.cat([protein_feat, glycan_feat], dim=1)
        x = self.fusion_fc(combined).unsqueeze(1)
        x = self.transformer_encoder(x)
        x = x.squeeze(1)
        return self.regressor(x)

# 训练函数
def train_model(model, train_loader, epochs=5):
    model.train()
    for epoch in range(epochs):
        total_loss = 0
        for protein, glycan, target in train_loader:
            protein = protein.to(device)
            glycan = glycan.to(device)
            target = target.to(device)
            
            optimizer.zero_grad()
            outputs = model(protein, glycan).squeeze()
            loss = criterion(outputs, target)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        print(f'Epoch {epoch+1}, Loss: {total_loss/len(train_loader):.4f}')

# 评估函数
def evaluate_model(model, test_loader):
    model.eval()
    predictions, truths = [], []
    with torch.no_grad():
        for protein, glycan, target in test_loader:
            protein = protein.to(device)
            glycan = glycan.to(device)
            outputs = model(protein, glycan).squeeze().cpu()
            predictions.extend(outputs.numpy())
            truths.extend(target.numpy())
    
    print(f'MAE: {mean_absolute_error(truths, predictions):.4f}')
    print(f'MSE: {mean_squared_error(truths, predictions):.4f}')
    print(f'R²: {r2_score(truths, predictions):.4f}')
    return predictions, truths

if __name__ == "__main__":
    # 初始化组件
    protein_tensor, glycan_tensor, target_tensor = load_data()
    dataset = CustomDataset(protein_tensor, glycan_tensor, target_tensor)
    
    # 数据集划分
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    train_dataset, test_dataset = random_split(dataset, [train_size, test_size])
    
    # 数据加载器
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True, pin_memory=True, num_workers=4)
    test_loader = DataLoader(test_dataset, batch_size=32)
    
    # 模型初始化
    protein_extractor = ProteinFeatureExtractor().to(device)
    glycan_extractor = GlycanFeatureExtractor().to(device)
    model = TransformerRegression(protein_extractor, glycan_extractor).to(device)
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
    
    # 训练与评估
    train_model(model, train_loader)
    preds, truths = evaluate_model(model, test_loader)
    
    # 残差图
    residuals = np.array(truths) - np.array(preds)
    plt.figure(figsize=(10,6))
    plt.scatter(preds, residuals, alpha=0.5)
    plt.axhline(y=0, color='r', linestyle='--')
    plt.xlabel('Predicted Values')
    plt.ylabel('Residuals')
    plt.title('Residual Plot')
    plt.show()