import numpy as np
import pandas as pd
import math
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from torch.nn import TransformerEncoder, TransformerEncoderLayer
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import matplotlib.pyplot as plt

# 数据准备
combined_data = pd.read_csv("combined_encode_esm_no-z-score.csv", header=None)

# 特征分割
protein_features = combined_data.iloc[:, :1280].values  # ESM模型提取的蛋白质特征向量
glycan_encoder = combined_data.iloc[:, 1280:1280+342].values  # 寡糖单糖、二糖、三糖坐标系编码向量
targets = combined_data.iloc[:, -1].values  # RFU值

# 第一次分割：分出临时数据（80%）和测试集（20%）
(X_protein_temp, X_protein_test, 
 X_glycan_temp, X_glycan_test,
 y_temp, y_test) = train_test_split(protein_features, glycan_encoder, targets, 
                                  test_size=0.2, random_state=42)

# 第二次分割：从临时数据中分出训练集（75%的temp）和验证集（25%的temp）
# 最终比例：训练60%（0.8*0.75=0.6），验证20%（0.8*0.25=0.2），测试20%
(X_protein_train, X_protein_val,
 X_glycan_train, X_glycan_val,
 y_train, y_val) = train_test_split(X_protein_temp, X_glycan_temp, y_temp,
                                   test_size=0.25, random_state=42)  # 0.25*0.8=0.2

# 标准化：训练集fit，测试集仅transform
scaler_protein = StandardScaler()
X_protein_train = scaler_protein.fit_transform(X_protein_train)
X_protein_val = scaler_protein.transform(X_protein_val)  
X_protein_test = scaler_protein.transform(X_protein_test)  

scaler_glycan = StandardScaler()
X_glycan_train = scaler_glycan.fit_transform(X_glycan_train)
X_glycan_val = scaler_glycan.transform(X_glycan_val)     
X_glycan_test = scaler_glycan.transform(X_glycan_test)

# 保留原始RFU值（不改变原有变量名）
y_train_raw = y_train.copy()
y_val_raw = y_val.copy()
y_test_raw = y_test.copy()

# 对目标值进行标准化（仅用训练集数据拟合）
targets_scaler = StandardScaler()
y_train = targets_scaler.fit_transform(y_train.reshape(-1, 1)).flatten()
y_val = targets_scaler.transform(y_val.reshape(-1, 1)).flatten()
y_test = targets_scaler.transform(y_test.reshape(-1, 1)).flatten()

# 自定义数据集
class CustomDataset(Dataset):
    def __init__(self, protein, glycan, targets_scaled, targets_raw):  
        self.protein = torch.FloatTensor(protein)
        self.glycan = torch.FloatTensor(glycan)
        self.targets_scaled = torch.FloatTensor(targets_scaled)
        self.targets_raw = torch.FloatTensor(targets_raw)  # 存储原始值
        
    def __len__(self):
        return len(self.targets_scaled)
    
    def __getitem__(self, idx):
        return self.protein[idx], self.glycan[idx], self.targets_scaled[idx], self.targets_raw[idx]

# 创建数据加载器
batch_size = 128
train_dataset = CustomDataset(X_protein_train, X_glycan_train, y_train, y_train_raw)
val_dataset = CustomDataset(X_protein_val, X_glycan_val, y_val, y_val_raw)
test_dataset = CustomDataset(X_protein_test, X_glycan_test, y_test, y_test_raw)

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)  # 打乱顺序
val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

# 位置编码类
class PositionalEncoding(nn.Module):
    def __init__(self, d_model: int, dropout: float = 0.1, max_len: int = 342):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)
        position = torch.arange(max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
        pe = torch.zeros(max_len, d_model)
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x + self.pe[:x.size(1)]
        return self.dropout(x)

# 寡糖特征提取
class GlycanEncoder(nn.Module):
    def __init__(self, d_model=64, nhead=8, num_layers=2, dim_feedforward=256, dropout=0.1):
        super().__init__()
        
         # 投影层：每个糖类型计数→d_model
        self.glycan_proj = nn.Linear(1, d_model)
        
        # 位置编码
        self.pos_encoder = PositionalEncoding(d_model)
        
        # Transformer编码器
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model, nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=dropout, batch_first=True
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers)

    def forward(self, x):
        # 添加通道维度并投影
        x = x.unsqueeze(-1)          # (B,342,1)
        x = self.glycan_proj(x)      # (B,342,d_model)
        
        # 位置编码
        x = self.pos_encoder(x)
        
        # Transformer处理
        x = self.encoder(x)          # (B,342,d_model)
        
        # 全局平均池化
        return x.mean(dim=1) 

# Transformer模型类
class GlycanProteinTransformer(nn.Module):
    def __init__(self):
        super().__init__()

        # 糖链特征处理
        self.glycan_encoder = GlycanEncoder()

        # 蛋白质特征处理
        self.protein_fc = nn.Sequential(
            nn.Linear(1280, 256),
            nn.GELU(),
            nn.Dropout(0.3),
            nn.Linear(256, 64),
            nn.LayerNorm(64)
        )
        
        # 联合回归层
        self.regressor = nn.Sequential(
            nn.Linear(128, 64),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(64, 1)
        )
        
    def forward(self, protein, glycan):
        # 糖链特征
        g_feat = self.glycan_encoder(glycan)  # (B, 64)
        
        # 蛋白特征
        p_feat = self.protein_fc(protein)     # (B, 64)
        
        # 联合预测
        combined = torch.cat([g_feat, p_feat], dim=1)  # (B, 128)
        return self.regressor(combined).squeeze()

# 训练准备
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = GlycanProteinTransformer().to(device)
criterion = nn.MSELoss()
optimizer = optim.AdamW(model.parameters(), lr=0.005, weight_decay=1e-4)
scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=5, factor=0.5)

# 训练循环
num_epochs = 10
best_val_loss = float('inf')

for epoch in range(num_epochs):
    model.train()
    train_loss = 0.0
    for protein, glycan, targets_scaled, _ in train_loader:
        protein, glycan, targets_scaled = protein.to(device), glycan.to(device), targets_scaled.to(device)
        optimizer.zero_grad()
        outputs = model(protein, glycan)
        loss = criterion(outputs, targets_scaled)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        train_loss += loss.item() * protein.size(0)
    train_loss /= len(train_loader.dataset)
    
    model.eval()
    val_loss = 0.0
    with torch.no_grad():
        for protein, glycan, targets_scaled, _ in val_loader:
            protein, glycan, targets_scaled = protein.to(device), glycan.to(device), targets_scaled.to(device)
            outputs = model(protein, glycan)
            val_loss += criterion(outputs, targets_scaled).item() * protein.size(0)
    val_loss /= len(val_loader.dataset)
    scheduler.step(val_loss)
    
    print(f'Epoch {epoch+1:03d} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}')
    
    # 保存最佳模型
    if val_loss < best_val_loss:
        best_val_loss = val_loss
        torch.save(model.state_dict(), 'best_transformer_model.pth')

# 测试评估
model.load_state_dict(torch.load('best_transformer_model.pth'))
model.eval()
test_preds_scaled = []  # 标准化后的预测值
test_targets_scaled = []  # 标准化后的真实值
test_targets_raw = []     # 原始真实值

with torch.no_grad():
    for protein, glycan, targets_scaled, targets_raw in test_loader:  
        protein, glycan = protein.to(device), glycan.to(device)
        preds_scaled = model(protein, glycan).cpu().numpy()
        
        test_preds_scaled.extend(preds_scaled)
        test_targets_scaled.extend(targets_scaled.numpy())
        test_targets_raw.extend(targets_raw.numpy())

# 标准化尺度下的指标
mae_scaled = mean_absolute_error(test_targets_scaled, test_preds_scaled)
mse_scaled = mean_squared_error(test_targets_scaled, test_preds_scaled)

# 原始尺度下的指标
test_preds_raw = targets_scaler.inverse_transform(np.array(test_preds_scaled).reshape(-1, 1)).flatten()
mae_raw = mean_absolute_error(test_targets_raw, test_preds_raw)
mse_raw = mean_squared_error(test_targets_raw, test_preds_raw)

print(f'\nStandardized Scale Metrics:')
print(f'MAE (scaled): {mae_scaled:.4f}')
print(f'MSE (scaled): {mse_scaled:.4f}')

print(f'\nRaw Scale Metrics:')
print(f'MAE (raw): {mae_raw:.4f}')
print(f'MSE (raw): {mse_raw:.4f}')
print(f'R² (raw) : {r2_score(test_targets_raw, test_preds_raw):.4f}')

# 可视化结果
plt.figure(figsize=(8, 6))
plt.scatter(test_targets_raw, test_preds_raw, alpha=0.5)
plt.plot([min(test_targets_raw), max(test_targets_raw)], 
         [min(test_targets_raw), max(test_targets_raw)], 'r--')
plt.xlabel('True RFU (Raw)')
plt.ylabel('Predicted RFU (Raw)')
plt.title('Regression Performance on Raw RFU Scale')
plt.show()