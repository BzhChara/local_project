"""
修改日志:
1.回归模块添加归一化层
2.交叉注意力+梯度裁剪
3.残差函数添加
4.学习率与优化器修改
5.位置编码
"""
import math
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
combined_data = pd.read_csv("combined_encode_esm_new.csv", header=None)

# 特征分割
protein_features = combined_data.iloc[:, :1280].values
glycan_encoder = combined_data.iloc[:, 1280:1280+342].values
targets = combined_data.iloc[:, -1].values

# 数据分割
(X_protein_train, X_protein_test, 
 X_glycan_train, X_glycan_test,
 y_train, y_test) = train_test_split(protein_features, glycan_encoder, targets, 
                                   test_size=0.2, random_state=42)

# 标准化：训练集fit，测试集仅transform
scaler_protein = StandardScaler()
X_protein_train = scaler_protein.fit_transform(X_protein_train)
X_protein_test = scaler_protein.transform(X_protein_test)  

scaler_glycan = StandardScaler()
X_glycan_train = scaler_glycan.fit_transform(X_glycan_train)
X_glycan_test = scaler_glycan.transform(X_glycan_test)

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
batch_size = 32
train_dataset = CustomDataset(X_protein_train, X_glycan_train, y_train)
test_dataset = CustomDataset(X_protein_test, X_glycan_test, y_test)

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)  # 打乱顺序
test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

# 位置编码
class PositionalEncoding(nn.Module):
    def __init__(self, d_model: int, max_len: int = 5000):
        super().__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, d_model, 2).float() *
            (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)  # 偶数维度用sin
        pe[:, 1::2] = torch.cos(position * div_term)   # 奇数维度用cos
        pe = pe.unsqueeze(0).transpose(0, 1)           # [max_len, 1, d_model]
        self.register_buffer('pe', pe)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x + self.pe[:x.size(0), :]  # 自动广播到batch维度
        return x

# 交叉注意力层
class CrossAttention(nn.Module):
    def __init__(self, d_model, nhead):
        super().__init__()
        self.attn = nn.MultiheadAttention(d_model, nhead, batch_first=True)
        self.norm = nn.LayerNorm(d_model)
        self.dropout = nn.Dropout(0.2)  
    
    def forward(self, query, key_value):
        # 保存原始输入
        identity = query
        # 计算注意力
        attn_output, _ = self.attn(query, key_value, key_value)
        # 残差连接 + 归一化
        output = self.norm(identity + self.dropout(attn_output))
        return output

# 残差模块
class ResidualBlock(nn.Module):
    def __init__(self, in_dim, out_dim):
        super().__init__()
        self.block = nn.Sequential(
            nn.Linear(in_dim, out_dim),
            nn.GELU(),
            nn.LayerNorm(out_dim),
            nn.Dropout(0.3)
        )
        self.shortcut = nn.Linear(in_dim, out_dim) if in_dim != out_dim else nn.Identity()

    def forward(self, x):
        return self.block(x) + self.shortcut(x)
    
# 定义Transformer回归模型
class TransformerRegressor(nn.Module):
    def __init__(self, protein_dim=1280, glycan_dim=342, d_model=512, nhead=8, hidden_dim=2048):
        super().__init__()
        
        # 蛋白质投影层
        self.protein_proj = nn.Sequential(
            nn.LayerNorm(protein_dim),
            nn.Linear(protein_dim, d_model*2),
            nn.GELU(),
            ResidualBlock(d_model*2, d_model*2),
            nn.Linear(d_model*2, d_model)
        )

        # 寡糖投影层
        self.glycan_proj = nn.Sequential(
            nn.LayerNorm(glycan_dim),
            nn.Linear(glycan_dim, d_model//2),
            nn.GELU(),
            ResidualBlock(d_model//2, d_model//2),
            nn.Linear(d_model//2, d_model)
        )
        self.pos_encoder = PositionalEncoding(d_model)
        
        # 寡糖特征处理模块
        self.glycan_encoded = nn.TransformerEncoder(
            encoder_layer=nn.TransformerEncoderLayer(
                d_model=d_model, 
                    nhead=nhead, 
                    dim_feedforward=hidden_dim, 
                    batch_first=True
            ),
            num_layers=2
        )
        
        # 交叉注意力模块（带残差）
        self.glycan_to_protein_attn = CrossAttention(d_model, nhead)  # 糖→蛋白
        self.protein_to_glycan_attn = CrossAttention(d_model, nhead)  # 蛋白→糖     
            
        # 回归模块
        self.regressor = nn.Sequential(
            ResidualBlock(d_model*2, d_model),
            ResidualBlock(d_model, d_model//2),
            nn.Linear(d_model//2, 1)
        )

    def forward(self, protein, glycan):
        # 蛋白质特征处理
        protein_proj = self.protein_proj(protein).unsqueeze(1)  # [batch, 1, d_model]
        
        # 寡糖特征处理
        glycan_proj = self.glycan_proj(glycan).unsqueeze(1)     # [batch, 1, d_model]
        glycan_proj = self.pos_encoder(glycan_proj)  # 添加位置编码
        glycan_features = self.glycan_encoded(glycan_proj)
        
        # 交叉注意力
        # 糖特征关注蛋白质
        glycan_aware = self.glycan_to_protein_attn(glycan_features, protein_proj)
        # 蛋白质关注糖特征
        protein_aware = self.protein_to_glycan_attn(protein_proj, glycan_features)
        
        # 拼接并回归
        combined = torch.cat([glycan_aware.squeeze(1), protein_aware.squeeze(1)], dim=1)
        return self.regressor(combined).squeeze()

# 模型初始化
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = TransformerRegressor().to(device)
criterion = nn.MSELoss()
optimizer = optim.AdamW(model.parameters(), lr=5e-5, weight_decay=0.001)
scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(optimizer, T_0=20)

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

        # 梯度裁剪
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=2.5) 

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

    # 更新学习率
    scheduler.step() 
    
    print(f"Epoch {epoch+1}/{num_epochs} | Train Loss: {train_loss:.4f} | Test Loss: {test_loss:.4f} | LR: {scheduler.get_last_lr()[0]:.2e}")

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