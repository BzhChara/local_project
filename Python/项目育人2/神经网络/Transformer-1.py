import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import matplotlib.pyplot as plt

# 数据准备
combined_data = pd.read_csv("combined_encode_esm_z-score.csv", header=None)

# 特征分割
protein_features = combined_data.iloc[:, :1280].values
glycan_features = combined_data.iloc[:, 1280:1280+342].values
targets = combined_data.iloc[:, -1].values

# 数据标准化（蛋白质和寡糖特征）
scaler_protein = StandardScaler()
protein_features = scaler_protein.fit_transform(protein_features)

scaler_glycan = StandardScaler()
glycan_features = scaler_glycan.fit_transform(glycan_features)

# 数据分割
(X_protein_train, X_protein_test, 
 X_glycan_train, X_glycan_test,
 y_train, y_test) = train_test_split(protein_features, glycan_features, targets, 
                                   test_size=0.2, random_state=42)

# 自定义数据集
class CustomDataset(Dataset):
    def __init__(self, protein, glycan, targets):
        self.protein = torch.FloatTensor(protein)
        self.glycan = torch.FloatTensor(glycan)
        self.targets = torch.FloatTensor(targets)
        
    def __len__(self):
        return len(self.targets)
    
    def __getitem__(self, idx):
        return self.protein[idx], self.glycan[idx], self.targets[idx]

# 创建数据加载器
batch_size = 128
train_dataset = CustomDataset(X_protein_train, X_glycan_train, y_train)
test_dataset = CustomDataset(X_protein_test, X_glycan_test, y_test)

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

# 定义Transformer回归模型
class TransformerRegressor(nn.Module):
    def __init__(self, protein_dim=1280, glycan_dim=342, d_model=64, nhead=8, num_layers=2, dim_feedforward=256, dropout=0.1):
        super().__init__()
        
        # 寡糖特征处理模块
        self.glycan_proj = nn.Linear(1, d_model)
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model, nhead=nhead, 
            dim_feedforward=dim_feedforward, 
            dropout=dropout, batch_first=True
        )
        self.glycan_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        # 蛋白质特征处理
        self.protein_fc = nn.Sequential(
            nn.Linear(protein_dim, 256),
            nn.GELU(),
            nn.Dropout(0.2),
            nn.Linear(256, 64),
        )
        
        # 回归模块
        self.regressor = nn.Sequential(
            nn.Linear(2 * d_model, 128),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(128, 1)
        )
    
    def forward(self, protein, glycan):
        # 处理寡糖特征
        glycan = glycan.unsqueeze(-1)  # [batch, 342, 1]
        glycan_proj = self.glycan_proj(glycan)  # [batch, 342, d_model]
        glycan_encoded = self.glycan_encoder(glycan_proj)
        glycan_pool = torch.mean(glycan_encoded, dim=1)  # 全局池化
        
        # 处理蛋白质
        protein = self.protein_fc(protein)

        # 拼接特征
        combined = torch.cat([protein, glycan_pool], dim=1)
        
        # 回归预测
        return self.regressor(combined).squeeze()

# 模型初始化
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = TransformerRegressor().to(device)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=1e-4, weight_decay=1e-5)

# 训练过程
num_epochs = 20
train_losses, test_losses = [], []

for epoch in range(num_epochs):
    # 训练阶段
    model.train()
    epoch_loss = 0
    for protein, glycan, targets in train_loader:
        protein, glycan, targets = protein.to(device), glycan.to(device), targets.to(device)
        
        optimizer.zero_grad()
        outputs = model(protein, glycan)
        loss = criterion(outputs, targets)
        loss.backward()
        optimizer.step()
        
        epoch_loss += loss.item() * protein.size(0)
    train_loss = epoch_loss / len(train_loader.dataset)
    train_losses.append(train_loss)
    
    # 验证阶段
    model.eval()
    test_loss = 0
    with torch.no_grad():
        for protein, glycan, targets in test_loader:
            protein, glycan, targets = protein.to(device), glycan.to(device), targets.to(device)
            outputs = model(protein, glycan)
            test_loss += criterion(outputs, targets).item() * protein.size(0)
    test_loss = test_loss / len(test_loader.dataset)
    test_losses.append(test_loss)
    
    print(f"Epoch {epoch+1}/{num_epochs} | Train Loss: {train_loss:.4f} | Test Loss: {test_loss:.4f}")

# 模型评估
model.eval()
all_preds = []
all_targets = []
with torch.no_grad():
    for protein, glycan, targets in test_loader:
        protein, glycan = protein.to(device), glycan.to(device)
        preds = model(protein, glycan).cpu().numpy()
        all_preds.extend(preds)
        all_targets.extend(targets.numpy())

# 计算指标
mae = mean_absolute_error(all_targets, all_preds)
mse = mean_squared_error(all_targets, all_preds)
r2 = r2_score(all_targets, all_preds)

print(f"\nEvaluation Results:")
print(f"MAE: {mae:.4f}")
print(f"MSE: {mse:.4f}")
print(f"R²: {r2:.4f}")

# 绘制残差图
plt.figure(figsize=(10, 6))
residuals = np.array(all_targets) - np.array(all_preds)
plt.scatter(all_preds, residuals, alpha=0.5)
plt.xlabel("Predicted Values")
plt.ylabel("Residuals")
plt.title("Residual Plot")
plt.axhline(y=0, color='r', linestyle='--')
plt.grid(True)
plt.savefig("output1.png")
plt.close()

# 训练过程可视化
plt.figure(figsize=(10, 6))
plt.plot(train_losses, label='Training Loss')
plt.plot(test_losses, label='Validation Loss')
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("Training Process")
plt.legend()
plt.grid(True)
plt.savefig("output2.png")
plt.close()