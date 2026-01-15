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
protein_features = combined_data.iloc[:, :1280].values  # ESM模型提取的蛋白质特征向量
glycan_encoder = combined_data.iloc[:, 1280:1280+342].values  # 寡糖单糖、二糖、三糖坐标系编码向量
targets = combined_data.iloc[:, -1].values  # Z-score后的RFU值

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
batch_size = 256
train_dataset = CustomDataset(X_protein_train, X_glycan_train, y_train)
test_dataset = CustomDataset(X_protein_test, X_glycan_test, y_test)

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)  # 打乱顺序
test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

# 交叉注意力层
class CrossAttention(nn.Module):
    def __init__(self, d_model, nhead):
        super().__init__()
        self.attn = nn.MultiheadAttention(d_model, nhead, batch_first=True)
        self.out_proj = nn.Linear(d_model, d_model)
        self.norm = nn.LayerNorm(d_model)
        self.dropout = nn.Dropout(0.2)
    
    def forward(self, query, key_value):
        # 保存原始输入
        identity = query
        # 计算注意力
        attn_output, _ = self.attn(query, key_value, key_value)
        attn_output = self.out_proj(attn_output)
        # 残差连接 + 归一化
        output = self.norm(identity + self.dropout(attn_output))
        return output

# 糖链特征结构化处理
class GlycanStructureEncoder(nn.Module):
    def __init__(self, mono_dim=22, di_dim=154, tri_dim=166, output_dim=512):  
        super().__init__()
        self.block_dim = output_dim // 4

       # 各层投影
        self.mono_encoder = nn.Sequential(
            nn.Linear(mono_dim, self.block_dim),
            nn.GELU(),
            nn.LayerNorm(self.block_dim),
            nn.Dropout(0.3)
        )
        self.di_encoder = nn.Sequential(
            nn.Linear(di_dim, self.block_dim),
            nn.GELU(),
            nn.LayerNorm(self.block_dim),
            nn.Dropout(0.3)
        )
        self.tri_encoder = nn.Sequential(
            nn.Linear(tri_dim, self.block_dim),
            nn.GELU(),
            nn.LayerNorm(self.block_dim),
            nn.Dropout(0.3)
        )
        self.adaptive_proj = nn.Sequential(
            nn.Linear(self.block_dim, self.block_dim),
            nn.GELU()
        )

        # 位置编码
        self.pos_embed = nn.Parameter(torch.randn(3, self.block_dim))
        
        # 交互层
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=self.block_dim, 
            nhead=8,
            dim_feedforward=256,
            dropout=0.2,
            batch_first=True
        )
        self.hier_interaction = nn.TransformerEncoder(encoder_layer, num_layers=2)

    def forward(self, glycan):
        # 特征提取
        mono = self.mono_encoder(glycan[:, :22]).unsqueeze(1)
        di = self.di_encoder(glycan[:, 22:176]).unsqueeze(1)
        tri = self.tri_encoder(glycan[:, 176:342]).unsqueeze(1)
        
        # 拼接并添加位置编码
        levels = torch.cat([mono, di, tri], dim=1)
        levels += self.pos_embed.unsqueeze(0)
        
        # 层级交互
        interacted = self.hier_interaction(levels) 
        # 通过线性投影生成第四块（智能填充）
        new_block = self.adaptive_proj(interacted.mean(dim=1)).unsqueeze(1)
        interacted = torch.cat([interacted, new_block], dim=1)
        return interacted  # [B,4,D]

# 定义Transformer回归模型
class TransformerRegressor(nn.Module):
    def __init__(self, protein_dim=1280, glycan_dim=342, d_model=512, nhead=8):
        super().__init__()
        # 蛋白质分块
        self.protein_proj = nn.Sequential(
            nn.Linear(1280, d_model),
            nn.Unflatten(1, (4, d_model//4))  
        )

        # 糖链分块
        self.glycan_proj = GlycanStructureEncoder(output_dim=d_model)

        # 交叉注意力层
        self.cross_attn_layers = nn.ModuleList([
            CrossAttention(d_model=d_model//4, nhead=nhead)  
            for _ in range(4)  
        ])

        # 回归模块
        self.regressor = nn.Sequential(
            nn.Linear(2*d_model, 512),  
            nn.GELU(),
            nn.Dropout(0.3),
            nn.Linear(512, 256),
            nn.GELU(),
            nn.Linear(256, 1)
        )

    def forward(self, protein, glycan):
       # 特征投影
        protein_proj = self.protein_proj(protein)  # [B,4,64]
        glycan_proj = self.glycan_proj(glycan)      # [B,4,64]
        
        # 交替交叉注意力
        for attn_layer in self.cross_attn_layers:
            protein_proj = attn_layer(protein_proj, glycan_proj)
            glycan_proj = attn_layer(glycan_proj, protein_proj)
        
        # 拼接特征
        fused = torch.cat([
            protein_proj.flatten(1), 
            glycan_proj.flatten(1)
        ], dim=1)

        return self.regressor(fused).squeeze(-1)
    
# 模型初始化
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = TransformerRegressor().to(device)
criterion = nn.MSELoss()

# 训练过程
num_epochs = 100
train_losses, test_losses = [], []

# 优化器设置
optimizer = optim.Adam(model.parameters(), lr=5e-4)

# 带热身的余弦退火调度器
scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(
    optimizer,
    T_0=10,      # 周期长度
    T_mult=2,     # 周期倍增因子
    eta_min=1e-6  # 最小学习率
)

for epoch in range(num_epochs):
    # 训练阶段
    model.train()
    epoch_loss = 0
    for protein, glycan, targets in train_loader:
        protein = protein.to(device)
        glycan = glycan.to(device)
        targets = targets.to(device)
        
        # 梯度清零
        optimizer.zero_grad()

        # 前向传播
        outputs = model(protein, glycan)
        loss = criterion(outputs, targets)

        # 反向传播与裁剪
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
        
        # 参数更新
        optimizer.step()
        
        epoch_loss += loss.item() * protein.size(0)
    
    train_loss = epoch_loss / len(train_loader.dataset)
    train_losses.append(train_loss)
    
    # 验证阶段
    model.eval()
    test_loss = 0
    with torch.no_grad():
        for protein, glycan, targets in test_loader:
            protein = protein.to(device)
            glycan = glycan.to(device)
            targets = targets.to(device)
            outputs = model(protein, glycan)
            test_loss += criterion(outputs, targets).item() * protein.size(0)
    test_loss = test_loss / len(test_loader.dataset)
    test_losses.append(test_loss)
    
    # 获取当前学习率
    current_lr = optimizer.param_groups[0]['lr']  
    # 学习率调整
    scheduler.step()

    print(f"Epoch {epoch+1} | Train Loss: {train_loss:.4f} | Test Loss: {test_loss:.4f} | LR: {current_lr:.2e}")

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

# 训练过程可视化
plt.figure(figsize=(10, 6))
plt.plot(train_losses, label='Training Loss')
plt.plot(test_losses, label='Validation Loss')
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("Training Process")
plt.legend()
plt.grid(True)
plt.savefig("output1.png")
plt.close()

""""
# 绘制残差图
plt.figure(figsize=(10, 6))
residuals = np.array(all_targets) - np.array(all_preds)
plt.scatter(all_preds, residuals, alpha=0.5)
plt.xlabel("Predicted Values")
plt.ylabel("Residuals")
plt.title("Residual Plot")
plt.axhline(y=0, color='r', linestyle='--')
plt.grid(True)
plt.savefig("output2.png")
plt.close()
"""

