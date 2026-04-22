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
    def __init__(self, protein_extractor, glycan_extractor, d_model=128, nhead=8, num_layers=4):
        super().__init__()
        self.protein_extractor = protein_extractor
        self.glycan_extractor = glycan_extractor
        
        self.fusion_fc = nn.Linear(64, d_model)  
        encoder_layer = nn.TransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=512, batch_first=True)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        self.regressor = nn.Sequential(
            nn.Linear(d_model, 64),
            nn.Dropout(0.3),
            nn.ReLU(),
            nn.Linear(64, 1)
        )

    def forward(self, protein, glycan):
        p_feat = self.protein_extractor(protein).unsqueeze(1)  # (b,1,64)
        g_feat = self.glycan_extractor(glycan).unsqueeze(1)    # (b,1,64)
        combined = torch.cat([p_feat, g_feat], dim=1)          # (b,2,64)
        x = self.fusion_fc(combined)                           # (b,2,d_model)
        x = self.transformer_encoder(x)
        return self.regressor(x[:, 0])

# 训练函数
def train_model(model, train_loader, val_loader, optimizer, criterion, epochs, patience, model_path='best_model.pth'):
    best_val_loss = float('inf')
    counter = 0
    
    for epoch in range(epochs):
        # 训练阶段
        model.train()
        total_train_loss = 0.0
        for protein, glycan, target in train_loader:
            protein = protein.to(device)
            glycan = glycan.to(device)
            target = target.to(device)
            
            optimizer.zero_grad()
            outputs = model(protein, glycan).squeeze()
            loss = criterion(outputs, target)
            loss.backward()
            optimizer.step()
            total_train_loss += loss.item()
        
        # 验证阶段
        model.eval()
        total_val_loss = 0.0
        with torch.no_grad():
            for protein, glycan, target in val_loader:
                protein = protein.to(device)
                glycan = glycan.to(device)
                target = target.to(device)
                outputs = model(protein, glycan).squeeze()
                loss = criterion(outputs, target)
                total_val_loss += loss.item()
        
        avg_train_loss = total_train_loss / len(train_loader)
        avg_val_loss = total_val_loss / len(val_loader)

        scheduler.step(avg_val_loss)  # 根据验证损失调整学习率

        print(f'Epoch {epoch+1}: Train Loss {avg_train_loss:.4f}, Val Loss {avg_val_loss:.4f}', flush=True)

        # 早停判断
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            counter = 0
            torch.save(model.state_dict(), model_path)
        else:
            counter += 1
            if counter >= patience:
                print(f'Early stopping at epoch {epoch+1}', flush=True)
                break
    
    # 加载最佳模型
    model.load_state_dict(torch.load(model_path))
    return model

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
    
    print(f'MAE: {mean_absolute_error(truths, predictions):.4f}', flush=True)
    print(f'MSE: {mean_squared_error(truths, predictions):.4f}', flush=True)
    print(f'R²: {r2_score(truths, predictions):.4f}', flush=True)
    return predictions, truths

if __name__ == "__main__":
    # 初始化组件
    protein_tensor, glycan_tensor, target_tensor = load_data()
    dataset = CustomDataset(protein_tensor, glycan_tensor, target_tensor)
    
    # 数据集划分
    train_size = int(0.7 * len(dataset))
    val_size = int(0.1 * len(dataset))
    test_size = len(dataset) - train_size - val_size
    train_dataset, val_dataset, test_dataset = random_split(dataset, [train_size, val_size, test_size])
    
    # 数据加载器
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True, num_workers=4)
    val_loader = DataLoader(val_dataset, batch_size=32)
    test_loader = DataLoader(test_dataset, batch_size=32)
    
    # 模型初始化
    protein_extractor = ProteinFeatureExtractor().to(device)
    glycan_extractor = GlycanFeatureExtractor().to(device)
    model = TransformerRegression(protein_extractor, glycan_extractor).to(device)
    criterion = nn.MSELoss()

    # 改进的优化器配置（增加权重衰减）
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    # 动态学习率调度器配置
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, 
        mode='min',         # 监控验证损失最小值
        factor=0.5,         # 学习率衰减系数
        patience=5,         # 允许连续5个epoch无改善
        min_lr=1e-6         # 设置学习率下限
    )
    
    # 训练与评估
    model = train_model(model, train_loader, val_loader, optimizer, criterion, epochs=200, patience=5)
    preds, truths = evaluate_model(model, test_loader)
    
    # 残差图
    residuals = np.array(truths) - np.array(preds)
    plt.figure(figsize=(10,6))
    plt.scatter(preds, residuals, alpha=0.5)
    plt.axhline(y=0, color='r', linestyle='--')
    plt.xlabel('Predicted Values')
    plt.ylabel('Residuals')
    plt.title('Residual Plot')
    plt.savefig('result1.png')
    # plt.show()