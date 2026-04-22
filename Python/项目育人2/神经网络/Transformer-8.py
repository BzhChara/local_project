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
combined_data = pd.read_csv("combined_encode_esm_embedding_z-score.csv", header=None)

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

# 位置编码类
class PositionalEncoding(nn.Module):
    def __init__(self, d_model: int, dropout: float = 0.1, max_len: int = 5000):
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
        
        # 定义各糖类型参数
        self.section_dims = [22, 154, 166]  # 单糖/二糖/三糖维度
        self.type_embeddings = nn.Embedding(3, d_model//2)  # 类型嵌入
        
        # 子段特征处理器
        self.sub_encoders = nn.ModuleList([
            nn.Sequential(
                nn.Linear(1, d_model//2),  # 特征嵌入
                nn.GELU(),
                nn.LayerNorm(d_model//2)
            ) for _ in range(3)
        ])
        
        # Transformer编码器
        encoder_layer = TransformerEncoderLayer(
            d_model=d_model, nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=dropout, batch_first=True
        )
        self.transformer = TransformerEncoder(encoder_layer, num_layers)
        
        # 动态位置编码
        self.pos_encoder = nn.Parameter(torch.randn(1, sum(self.section_dims), d_model))

    def forward(self, glycan):
        # 拆分特征段
        sections = torch.split(glycan, self.section_dims, dim=1)  # 3个tensor
        
        # 分类型处理
        encoded = []
        for i, (section, encoder) in enumerate(zip(sections, self.sub_encoders)):
            # 子段特征嵌入 (B, N, 1) -> (B, N, d_model//2)
            feat = encoder(section.unsqueeze(-1)) 
            
            # 添加类型嵌入
            type_emb = self.type_embeddings(torch.tensor(i).to(glycan.device))  # (d_model)
            type_emb = type_emb.view(1, 1, -1).expand(feat.size(0), feat.size(1), -1)
            
            # 拼接特征和类型信息
            combined = torch.cat([feat, type_emb], dim=-1)  # (B, N, d_model)
            encoded.append(combined)
        
        # 合并所有糖单元
        full_sequence = torch.cat(encoded, dim=1)  # (B, 342, d_model)
        
        # 添加动态位置编码
        full_sequence += self.pos_encoder[:, :full_sequence.size(1), :]
        
        # Transformer编码
        return self.transformer(full_sequence)  # (B, 342, d_model)

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
        g_feat = self.glycan_encoder(glycan)
        g_pooled = g_feat.mean(dim=1)  # (B, d_model)
        
        # 蛋白特征
        p_feat = self.protein_fc(protein)    # (B, 64)
        
        # 联合预测
        combined = torch.cat([g_pooled, p_feat], dim=1)
        return self.regressor(combined).squeeze()

# 训练准备
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = GlycanProteinTransformer().to(device)
criterion = nn.MSELoss()
optimizer = optim.AdamW(model.parameters(), lr=0.005, weight_decay=1e-4)
scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=5, factor=0.5)

# 训练循环
num_epochs = 5
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
        torch.save(model.state_dict(), 'best_transformer_model.pth')

# 测试评估
model.load_state_dict(torch.load('best_transformer_model.pth'))
model.eval()
test_preds, test_targets = [], []

with torch.no_grad():
    for protein, glycan, targets in test_loader:
        protein, glycan = protein.to(device), glycan.to(device)
        preds = model(protein, glycan).cpu().numpy()
        test_preds.extend(preds)
        test_targets.extend(targets.numpy())

# 计算指标
mae = mean_absolute_error(test_targets, test_preds)
mse = mean_squared_error(test_targets, test_preds)
r2 = r2_score(test_targets, test_preds)

print(f'\nTest Performance:')
print(f'MAE: {mae:.4f}')
print(f'MSE: {mse:.4f}')
print(f'R²:  {r2:.4f}')

# 可视化结果
plt.figure(figsize=(8, 6))
plt.scatter(test_targets, test_preds, alpha=0.5)
plt.plot([min(test_targets), max(test_targets)], [min(test_targets), max(test_targets)], 'r--')
plt.xlabel('True Z-score')
plt.ylabel('Predicted Z-score')
plt.title('Regression Performance on Test Set')
plt.show()