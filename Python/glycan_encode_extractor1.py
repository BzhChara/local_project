import torch
import torch.nn as nn
import math

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

class GlycanEncoder(nn.Module):
    def __init__(self, d_model=64, nhead=8, num_layers=2, dim_feedforward=256, dropout=0.1):
        super().__init__()
        self.glycan_proj = nn.Linear(1, d_model)  # 输入为342维，每个位置编码为1维
        self.pos_encoder = PositionalEncoding(d_model)
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model, nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=dropout, batch_first=True
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers)
        
    def forward(self, x):
        x = x.unsqueeze(-1)          # (B, 342, 1)
        x = self.glycan_proj(x)      # (B, 342, d_model)
        x = self.pos_encoder(x)
        x = self.encoder(x)          # (B, 342, d_model)
        return x.mean(dim=1)         # 全局平均池化 → (B, d_model)
    
class JointModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.glycan_encoder = GlycanEncoder()  # 输出64维
        
        # 蛋白质特征处理层
        self.protein_fc = nn.Sequential(
            nn.Linear(1280, 256),
            nn.GELU(),
            nn.Dropout(0.3),
            nn.Linear(256, 64),
            nn.LayerNorm(64)
        )
        
        # 联合回归层
        self.regressor = nn.Sequential(
            nn.Linear(128, 64),     # 64（糖） + 64（蛋白）= 128
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(64, 1)
        )
        
    def forward(self, protein, glycan):
        g_feat = self.glycan_encoder(glycan)   # (B, 64)
        p_feat = self.protein_fc(protein)      # (B, 64)
        combined = torch.cat([g_feat, p_feat], dim=1)
        return self.regressor(combined).squeeze()
    
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from torch.utils.data import Dataset, DataLoader

# 加载数据
combined_data = pd.read_csv("combined_encode_esm_no-z-score.csv", header=None)
protein_features = combined_data.iloc[:, :1280].values
glycan_encoder = combined_data.iloc[:, 1280:1280+342].values
targets = combined_data.iloc[:, -1].values

# 分割数据集（训练:验证:测试 = 60%:20%:20%）
X_protein_train, X_protein_test, X_glycan_train, X_glycan_test, y_train, y_test = train_test_split(
    protein_features, glycan_encoder, targets, test_size=0.2, random_state=42
)
X_protein_train, X_protein_val, X_glycan_train, X_glycan_val, y_train, y_val = train_test_split(
    X_protein_train, X_glycan_train, y_train, test_size=0.25, random_state=42
)

# 标准化处理（注意：仅用训练集拟合）
scaler_protein = StandardScaler()
X_protein_train = scaler_protein.fit_transform(X_protein_train)
X_protein_val = scaler_protein.transform(X_protein_val)
X_protein_test = scaler_protein.transform(X_protein_test)

scaler_glycan = StandardScaler()
X_glycan_train = scaler_glycan.fit_transform(X_glycan_train)
X_glycan_val = scaler_glycan.transform(X_glycan_val)
X_glycan_test = scaler_glycan.transform(X_glycan_test)

# 目标值标准化
target_scaler = StandardScaler()
y_train = target_scaler.fit_transform(y_train.reshape(-1, 1)).flatten()
y_val = target_scaler.transform(y_val.reshape(-1, 1)).flatten()
y_test = target_scaler.transform(y_test.reshape(-1, 1)).flatten()

# 自定义数据集类
class GlycanDataset(Dataset):
    def __init__(self, protein, glycan, targets):
        self.protein = torch.FloatTensor(protein)
        self.glycan = torch.FloatTensor(glycan)
        self.targets = torch.FloatTensor(targets)
    
    def __len__(self):
        return len(self.targets)
    
    def __getitem__(self, idx):
        return self.protein[idx], self.glycan[idx], self.targets[idx]

# 创建DataLoader
batch_size = 64
train_dataset = GlycanDataset(X_protein_train, X_glycan_train, y_train)
val_dataset = GlycanDataset(X_protein_val, X_glycan_val, y_val)
test_dataset = GlycanDataset(X_protein_test, X_glycan_test, y_test)

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = JointModel().to(device)
criterion = nn.MSELoss()
optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=5, factor=0.5)

# 训练循环
num_epochs = 50
best_val_loss = float('inf')

for epoch in range(num_epochs):
    model.train()
    train_loss = 0.0
    for protein, glycan, targets in train_loader:
        protein, glycan, targets = protein.to(device), glycan.to(device), targets.to(device)
        optimizer.zero_grad()
        outputs = model(protein, glycan)
        loss = criterion(outputs, targets)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        train_loss += loss.item() * protein.size(0)
    train_loss /= len(train_loader.dataset)
    
    model.eval()
    val_loss = 0.0
    with torch.no_grad():
        for protein, glycan, targets in val_loader:
            protein, glycan, targets = protein.to(device), glycan.to(device), targets.to(device)
            outputs = model(protein, glycan)
            val_loss += criterion(outputs, targets).item() * protein.size(0)
    val_loss /= len(val_loader.dataset)
    scheduler.step(val_loss)
    
    print(f'Epoch {epoch+1:03d} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}')
    
    # 保存最佳模型
    if val_loss < best_val_loss:
        best_val_loss = val_loss
        torch.save(model.glycan_encoder.state_dict(), 'best_glycan_encoder.pth')  # 仅保存GlycanEncoder

# 加载训练好的GlycanEncoder
glycan_encoder = GlycanEncoder().to(device)
glycan_encoder.load_state_dict(torch.load('best_glycan_encoder.pth'))
glycan_encoder.eval()

# 定义特征提取函数
def extract_features(glycan_data, glycan_scaler):
    glycan_scaled = glycan_scaler.transform(glycan_data)
    glycan_tensor = torch.FloatTensor(glycan_scaled).to(device)
    with torch.no_grad():
        features = glycan_encoder(glycan_tensor).cpu().numpy()
    return features

# 提取各数据集的特征
glycan_train_feats = extract_features(X_glycan_train, scaler_glycan)
glycan_val_feats = extract_features(X_glycan_val, scaler_glycan)
glycan_test_feats = extract_features(X_glycan_test, scaler_glycan)

# 合并蛋白质与寡糖特征
def merge_save_features(protein_feats, glycan_feats, targets, filename):
    merged = np.concatenate([protein_feats, glycan_feats], axis=1)
    merged_df = pd.DataFrame(merged)
    merged_df['RFU'] = targets
    merged_df.to_csv(filename, index=False, header=False)

merge_save_features(X_protein_train, glycan_train_feats, y_train, 'train_h2o.csv')
merge_save_features(X_protein_val, glycan_val_feats, y_val, 'val_h2o.csv')
merge_save_features(X_protein_test, glycan_test_feats, y_test, 'test_h2o.csv')
