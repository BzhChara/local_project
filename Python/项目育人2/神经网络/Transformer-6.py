import numpy as np
import pandas as pd
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
combined_data = pd.read_csv("merged_features.csv", header=None)

# 特征分割
protein_features = combined_data.iloc[:, :1280].values  # ESM模型提取的蛋白质特征向量
glycan_encoder = combined_data.iloc[:, 1280:1280+342].values  # 寡糖单糖、二糖、三糖坐标系编码向量
targets = combined_data.iloc[:, -1].values  # Z-score后的RFU值

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
val_dataset = CustomDataset(X_protein_val, X_glycan_val, y_val)
test_dataset = CustomDataset(X_protein_test, X_glycan_test, y_test)

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)  # 打乱顺序
val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=5000):
        super().__init__()
        self.dropout = nn.Dropout(p=0.1)
        
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-torch.log(torch.tensor(10000.0)) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        self.register_buffer('pe', pe)

    def forward(self, x):
        x = x + self.pe[:, :x.size(1)]
        return self.dropout(x)

class ProteinEncoder(nn.Module):
    def __init__(self):
        super().__init__()
        # 分块处理(1280 -> 80x16)
        self.chunk_embedding = nn.Linear(16, 64)
        self.pos_encoder = PositionalEncoding(64)
        encoder_layer = TransformerEncoderLayer(d_model=64, nhead=8, dim_feedforward=256, dropout=0.1, batch_first=True)
        self.transformer_encoder = TransformerEncoder(encoder_layer, num_layers=3)
        self.pool = nn.AdaptiveAvgPool1d(1)

    def forward(self, x):
        # x shape: (batch_size, 1280)
        x = x.view(x.size(0), 80, 16)  # 分块为80个16维特征
        x = self.chunk_embedding(x)    # (batch, 80, 64)
        x = self.pos_encoder(x)
        x = self.transformer_encoder(x)  # (batch, 80, 64)
        x = x.transpose(1, 2)  # (batch, 64, 80)
        x = self.pool(x).squeeze(2)  # (batch, 64)
        return x
        
class GlycanEncoder(nn.Module):
    def __init__(self):
        super().__init__()
        # 分层特征嵌入
        self.mono_encoder = nn.Sequential(
            nn.Linear(22, 64),
            nn.LayerNorm(64),
            nn.GELU()
        )
        
        self.di_encoder = nn.Sequential(
            nn.Linear(154, 128),
            nn.LayerNorm(128),
            nn.GELU(),
            nn.Linear(128, 64),
            nn.Dropout(0.2)
        )
        
        self.tri_encoder = nn.Sequential(
            nn.Linear(166, 128),
            nn.LayerNorm(128),
            nn.GELU(),
            nn.Linear(128, 64),
            nn.Dropout(0.2)
        )
        
        # 注意力融合层
        self.attention = nn.MultiheadAttention(64, 4, batch_first=True)
        self.pool = nn.AdaptiveAvgPool1d(1)

    def forward(self, x):
        # 分解特征
        mono = x[:, :22]       # (bs,22)
        di = x[:, 22:176]      # (bs,154) 
        tri = x[:, 176:342]    # (bs,166)
        
        # 分层编码
        mono_feat = self.mono_encoder(mono).unsqueeze(1)  # (bs,1,64)
        di_feat = self.di_encoder(di).unsqueeze(1)        # (bs,1,64)
        tri_feat = self.tri_encoder(tri).unsqueeze(1)     # (bs,1,64)
        
        # 拼接特征序列
        features = torch.cat([mono_feat, di_feat, tri_feat], dim=1)  # (bs,3,64)
        
        # 注意力融合
        attn_out, _ = self.attention(features, features, features)  # (bs,3,64)
        
        # 全局特征提取
        return self.pool(attn_out.transpose(1,2)).squeeze(2)  # (bs,64)

# 定义Transformer回归模型
class TransformerRegressor(nn.Module):
    def __init__(self):
        super().__init__()
        self.protein_encoder = ProteinEncoder()
        self.glycan_encoder = GlycanEncoder()
        
        # 交叉注意力机制
        self.cross_attn = nn.MultiheadAttention(embed_dim=64, num_heads=4, batch_first=True)

        # 动态特征融合
        self.regressor = nn.Sequential(
            nn.Linear(192, 128),
            nn.BatchNorm1d(128),
            nn.GELU(),
            nn.Dropout(0.3),
            nn.Linear(128, 64),
            nn.LayerNorm(64),
            nn.GELU(),
            nn.Dropout(0.2),
            nn.Linear(64, 1)
        )

    def forward(self, protein, glycan):
        prot_feat = self.protein_encoder(protein)
        gly_feat = self.glycan_encoder(glycan)

        # 交叉注意力
        cross_feat, _ = self.cross_attn(
            query=prot_feat.unsqueeze(1),
            key=gly_feat.unsqueeze(1),
            value=gly_feat.unsqueeze(1))

        combined = torch.cat([
            prot_feat, 
            gly_feat,
            cross_feat.squeeze(1)], dim=1)  # (bs,128+128+128=384)
            
        return self.regressor(combined)
    
# 模型训练配置
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = TransformerRegressor().to(device)
criterion = nn.MSELoss() 
optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4, weight_decay=1e-5)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.5, patience=5)

# 早停机制
class EarlyStopper:
    def __init__(self, patience=10, min_delta=0):
        self.patience = patience
        self.min_delta = min_delta
        self.counter = 0
        self.min_validation_loss = float('inf')

    def __call__(self, validation_loss):
        if validation_loss < self.min_validation_loss - self.min_delta:
            self.min_validation_loss = validation_loss
            self.counter = 0
        else:
            self.counter += 1
            if self.counter >= self.patience:
                return True
        return False

early_stopper = EarlyStopper(patience=15, min_delta=0.001)

# 训练函数
def train_epoch(model, loader):
    model.train()
    total_loss = 0
    for protein, glycan, targets in loader:
        protein = protein.to(device)
        glycan = glycan.to(device)
        targets = targets.to(device).unsqueeze(1)
        
        optimizer.zero_grad()
        outputs = model(protein, glycan)
        loss = criterion(outputs, targets)
        loss.backward()
        
        # 梯度裁剪
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        
        optimizer.step()
        total_loss += loss.item() * protein.size(0)
    return total_loss / len(loader.dataset)

# 验证函数
def validate(model, loader):
    model.eval()
    total_loss = 0
    with torch.no_grad():
        for protein, glycan, targets in loader:
            protein = protein.to(device)
            glycan = glycan.to(device)
            targets = targets.to(device).unsqueeze(1)
            
            outputs = model(protein, glycan)
            loss = criterion(outputs, targets)
            total_loss += loss.item() * protein.size(0)
    return total_loss / len(loader.dataset)

# 测试函数
def test(model, loader):
    model.eval()
    predictions = []
    truths = []
    with torch.no_grad():
        for protein, glycan, targets in loader:
            protein = protein.to(device)
            glycan = glycan.to(device)
            
            outputs = model(protein, glycan)
            predictions.append(outputs.cpu().numpy())
            truths.append(targets.numpy())
    
    y_pred = np.concatenate(predictions)
    y_true = np.concatenate(truths)
    
    metrics = {
        'MSE': mean_squared_error(y_true, y_pred),
        'MAE': mean_absolute_error(y_true, y_pred),
        'R2': r2_score(y_true, y_pred)
    }
    return metrics

# 主训练循环
best_val_loss = float('inf')
for epoch in range(50):
    train_loss = train_epoch(model, train_loader)
    val_loss = validate(model, val_loader)
    
    # 学习率调整
    scheduler.step(val_loss)
    
    # 保存最佳模型
    if val_loss < best_val_loss:
        best_val_loss = val_loss
        torch.save(model.state_dict(), 'best_model.pth')
    
    # 打印训练信息
    print(f'Epoch {epoch:03d} | '
          f'Train Loss: {train_loss:.4f} | '
          f'Val Loss: {val_loss:.4f} | '
          f'LR: {optimizer.param_groups[0]["lr"]:.2e}')
    
    # 早停检查
    if early_stopper(val_loss):
        print("Early stopping triggered!")
        break

# 最终测试
model.load_state_dict(torch.load('best_model.pth', weights_only=True))
test_metrics = test(model, test_loader)
print("\nTest Results:")
print(f"MSE: {test_metrics['MSE']:.4f}")
print(f"MAE: {test_metrics['MAE']:.4f}")
print(f"R²: {test_metrics['R2']:.4f}")