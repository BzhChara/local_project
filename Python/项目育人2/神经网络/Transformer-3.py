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
batch_size = 64
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

# 层次化交叉注意力设计
class HierarchicalCrossAttention(nn.Module):
    def __init__(self, d_model, nhead=8):
        super().__init__()
        # 蛋白质特征投影
        self.protein_proj = nn.Sequential(ResidualBlock(d_model//3, d_model//3))

        # 三层独立注意力机制
        self.mono_attn = CrossAttention(d_model//3, nhead)
        self.di_attn = CrossAttention(d_model//3, nhead)
        self.tri_attn = CrossAttention(d_model//3, nhead)

        # 融合层
        self.fusion = nn.Sequential(
            ResidualBlock(128, 256),
            ResidualBlock(256, 128),
            nn.Dropout(0.2)
        )
        
        # 最终自注意力层
        self.final_attn = nn.MultiheadAttention(
            embed_dim=d_model//3,  
            num_heads=8,
            batch_first=True
        )
        self.final_norm = nn.LayerNorm(d_model//3)

    def forward(self, protein_chunks, glycan_hierarchy):
        # 投影蛋白质特征
        protein_proj = self.protein_proj(protein_chunks)  # [B,31,128]

        # 输入维度验证
        # print(f"Protein feats: {protein_proj.shape}")  # 应为[batch, 31, 128]
        # print(f"Glycan feats: {glycan_hierarchy.shape}")    # 应为[batch, 3, 128]
        
        # 拆分糖层级特征
        mono, di, tri = torch.chunk(glycan_hierarchy, 3, dim=1)
        
        # 各层级与蛋白质交互
        mono_out = self.mono_attn(mono, protein_proj)
        di_out = self.di_attn(di, protein_proj)
        tri_out = self.tri_attn(tri, protein_proj)
        
        # 拼接并融合
        combined = torch.cat([mono_out, di_out, tri_out], dim=1)
        
        # 最终自注意力
        attn_output, _ = self.final_attn(combined, combined, combined)
        output = self.final_norm(combined + attn_output)
        return self.fusion(output)
    
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
    
# 蛋白质特征分块策略
class ProteinChunkEncoder(nn.Module):
    def __init__(self, protein_dim=1280, chunk_size=80, d_model=384):
        super().__init__()
        self.num_chunks = (protein_dim - chunk_size) // (chunk_size//2) + 1  # 重叠分块计算
        self.chunk_proj = nn.Sequential(
            nn.Linear(chunk_size, d_model//3),
            nn.GELU(),
            nn.LayerNorm(d_model//3),
            nn.Dropout(0.2)
        )
        
    def forward(self, protein):
        # 使用滑动窗口分块 (50%重叠)
        chunks = protein.unfold(1, 80, 40)  # [batch, num_chunks, 80]
        
        # 投影到统一维度
        return self.chunk_proj(chunks)  # [batch, num_chunks, d_model]

# 糖链特征结构化处理
class GlycanStructureEncoder(nn.Module):
    def __init__(self, mono_dim=22, di_dim=154, tri_dim=166, d_model=384):
        super().__init__()
        # 分层投影
        self.mono_encoder = nn.Sequential(
            nn.Linear(mono_dim, d_model//3),
            nn.GELU(),
            nn.LayerNorm(d_model//3),
            nn.Dropout(0.2)
        )
        
        self.di_encoder = nn.Sequential(
            nn.Linear(di_dim, d_model//3),
            nn.GELU(),
            nn.LayerNorm(d_model//3),
            nn.Dropout(0.2)
        )
        
        self.tri_encoder = nn.Sequential(
            nn.Linear(tri_dim, d_model//3),
            nn.GELU(),
            nn.LayerNorm(d_model//3),
            nn.Dropout(0.2)
        )
        
        # 层级位置编码 + 类型嵌入
        self.mono_pos = nn.Parameter(torch.randn(1, 1, d_model//3))
        self.di_pos = nn.Parameter(torch.randn(1, 1, d_model//3))
        self.tri_pos = nn.Parameter(torch.randn(1, 1, d_model//3))

         # 深层交互模块
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model//3,
            nhead=8,
            dim_feedforward=256,
            dropout=0.2,
            batch_first=True
        )
        self.hier_interaction = nn.TransformerEncoder(
            encoder_layer,
            num_layers=2 
        )
        self.attn_norm = nn.LayerNorm(d_model//3)

    def forward(self, glycan):
        # 特征提取与增强
        mono = self.mono_encoder(glycan[:, :22]).unsqueeze(1) 
        mono = mono + self.mono_pos
        
        di = self.di_encoder(glycan[:, 22:176]).unsqueeze(1)
        di = di + self.di_pos
        
        tri = self.tri_encoder(glycan[:, 176:342]).unsqueeze(1)
        tri = tri + self.tri_pos

        # 拼接层级特征 [batch, 3, d_model//3]
        levels = torch.cat([mono, di, tri], dim=1)

        # 深层交互处理（包含自注意力与前馈网络）
        interacted = self.hier_interaction(levels)
        
        # 残差连接+归一化
        return self.attn_norm(levels + interacted)
    
# 定义Transformer回归模型
class TransformerRegressor(nn.Module):
    def __init__(self, protein_dim=1280, glycan_dim=342, d_model=384, nhead=8):
        super().__init__()
        
        # 蛋白质处理
        self.protein_encoder = ProteinChunkEncoder(protein_dim, d_model=d_model)

        # 糖链处理
        self.glycan_encoder = GlycanStructureEncoder(d_model=d_model)

        # 交叉注意力
        self.hier_cross_attn = HierarchicalCrossAttention(d_model, nhead)
            
        # 回归模块
        self.regressor = nn.Sequential(
            ResidualBlock(128, 512),  
            ResidualBlock(512, 256),
            ResidualBlock(256, 128),
            ResidualBlock(128, 64),
            nn.Linear(64, 1)
        )

    def forward(self, protein, glycan):
        # 打印输入维度
        # print(f"Input protein shape: {protein.shape}")  # 应为[batch, 1280]
        # print(f"Input glycan shape: {glycan.shape}")   # 应为[batch, 342]

        # 特征编码
        protein_feats = self.protein_encoder(protein)  # [B, C, D]
        glycan_feats = self.glycan_encoder(glycan)     # [B, 3, D/3]
        
        # 打印编码后维度
        # print(f"Encoded protein: {protein_feats.shape}")  # 应为[batch, 31, 128]
        # print(f"Encoded glycan: {glycan_feats.shape}")    # 应为[batch, 3, 128]

        # 层次化交叉注意力
        combined = self.hier_cross_attn(protein_feats, glycan_feats)
        # print(f"Combined shape: {combined.shape}")  # 应为[batch, 3, 128]
        
        # 全局平均后回归
        return self.regressor(combined.mean(dim=1)).squeeze(-1)

# 模型初始化
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = TransformerRegressor().to(device)
criterion = nn.MSELoss()

# 训练过程
num_epochs = 50
train_losses, test_losses = [], []

# 优化器设置
optimizer = optim.AdamW(model.parameters(), lr=1e-4, weight_decay=0.01)

# 学习率调度器配置
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
    optimizer, 
    mode='min',      # 监控验证损失的最小值
    factor=0.5,      # 学习率衰减因子
    patience=3,      # 容忍3个epoch没有改进
    threshold=0.01, # 超过此阈值才视为改进
    cooldown=1       # 调整后冷却1个epoch
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
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=2.5)
        
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
    
    # 学习率调整
    scheduler.step(test_loss)

    current_lr = optimizer.param_groups[0]['lr']  # 获取当前学习率
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

