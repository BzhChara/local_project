import torch
import torch.nn as nn
import torch.optim as optim
import math
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from torch.utils.data import Dataset, DataLoader
import matplotlib.pyplot as plt

# ==================== 数据准备 ====================
combined_data = pd.read_csv("combined_encode_esm_new.csv", header=None)

# 特征分割
protein_features = combined_data.iloc[:, :1280].values
glycan_encoder = combined_data.iloc[:, 1280:1280+342].values
targets = combined_data.iloc[:, -1].values.reshape(-1, 1)

# 数据检查
print(f"蛋白质特征形状: {protein_features.shape}")
print(f"寡糖特征形状: {glycan_encoder.shape}")
print(f"目标值形状: {targets.shape}")

# 检查NaN值
assert not np.isnan(protein_features).any(), "蛋白质特征包含NaN值"
assert not np.isnan(glycan_encoder).any(), "寡糖特征包含NaN值"
assert not np.isnan(targets).any(), "目标值包含NaN值"

# 全局索引打乱
indices = np.arange(len(targets))
np.random.seed(42)
shuffled_indices = np.random.permutation(indices)

# 按打乱后的索引重新排列数据
protein_shuffled = protein_features[shuffled_indices]
glycan_shuffled = glycan_encoder[shuffled_indices]
targets_shuffled = targets[shuffled_indices]

# 划分训练/验证/测试 (70-15-15)
X_protein_train, X_protein_temp, X_glycan_train, X_glycan_temp, y_train, y_temp = train_test_split(
    protein_shuffled, glycan_shuffled, targets_shuffled, 
    test_size=0.3, random_state=42
)

X_protein_val, X_protein_test, X_glycan_val, X_glycan_test, y_val, y_test = train_test_split(
    X_protein_temp, X_glycan_temp, y_temp,
    test_size=0.5, random_state=42
)

# 标准化
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

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

# ==================== 模型定义 ====================
class MultiHeadAttention(nn.Module):
    def __init__(self, d_model, num_heads):
        super().__init__()
        assert d_model % num_heads == 0, "d_model必须能被num_heads整除"
        
        self.d_model = d_model
        self.num_heads = num_heads
        self.d_k = d_model // num_heads
        
        self.W_q = nn.Linear(d_model, d_model)
        self.W_k = nn.Linear(d_model, d_model)
        self.W_v = nn.Linear(d_model, d_model)
        self.W_o = nn.Linear(d_model, d_model)
        
    def scaled_dot_product_attention(self, Q, K, V):
        attn_scores = torch.matmul(Q, K.transpose(-2, -1)) / math.sqrt(self.d_k)
        attn_probs = torch.softmax(attn_scores, dim=-1)
        return torch.matmul(attn_probs, V)
        
    def split_heads(self, x):
        return x.view(x.size(0), x.size(1), self.num_heads, self.d_k).transpose(1, 2)
        
    def combine_heads(self, x):
        return x.transpose(1, 2).contiguous().view(x.size(0), x.size(2), self.d_model)
        
    def forward(self, x):
        Q = self.split_heads(self.W_q(x))
        K = self.split_heads(self.W_k(x))
        V = self.split_heads(self.W_v(x))
        attn_output = self.scaled_dot_product_attention(Q, K, V)
        return self.W_o(self.combine_heads(attn_output))

class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_seq_len=342):
        super().__init__()
        position = torch.arange(max_seq_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
        
        pe = torch.zeros(max_seq_len, d_model)
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe.unsqueeze(0))
        
    def forward(self, x):
        return x + self.pe[:, :x.size(1)]

class PositionWiseFFN(nn.Module):
    def __init__(self, d_model, d_ff):
        super().__init__()
        self.fc1 = nn.Linear(d_model, d_ff)
        self.fc2 = nn.Linear(d_ff, d_model)
        self.relu = nn.ReLU()
        
    def forward(self, x):
        return self.fc2(self.relu(self.fc1(x)))

class TransformerRegressor(nn.Module):
    def __init__(self, protein_dim=1280, glycan_dim=342, 
                 d_model=512, num_heads=8, num_layers=4, 
                 d_ff=2048, dropout=0.1):
        super().__init__()
        
        # 蛋白质特征处理
        self.protein_net = nn.Sequential(
            nn.Linear(protein_dim, d_model),
            nn.LayerNorm(d_model),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # 寡糖特征处理
        self.glycan_embed = nn.Linear(1, d_model)  # 每个特征视为一个token
        self.pos_encoder = PositionalEncoding(d_model)
        self.attn_layers = nn.ModuleList([
            MultiHeadAttention(d_model, num_heads) for _ in range(num_layers)
        ])
        self.ffn = PositionWiseFFN(d_model, d_ff)
        self.layer_norm = nn.LayerNorm(d_model)
        
        # 回归层
        self.regressor = nn.Sequential(
            nn.Linear(2*d_model, d_model),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(d_model, 1)
        )

    def forward(self, protein, glycan):
        # 蛋白质特征
        protein_feat = self.protein_net(protein)
        
        # 寡糖特征处理
        glycan = glycan.unsqueeze(-1)  # [batch, 342, 1]
        glycan_embed = self.glycan_embed(glycan)  # [batch, 342, d_model]
        glycan_embed = self.pos_encoder(glycan_embed)
        
        # 自注意力处理
        for attn_layer in self.attn_layers:
            glycan_embed = attn_layer(glycan_embed) + glycan_embed
            glycan_embed = self.ffn(glycan_embed) + glycan_embed
            glycan_embed = self.layer_norm(glycan_embed)
        
        # 全局平均池化
        glycan_feat = torch.mean(glycan_embed, dim=1)
        
        # 特征融合
        combined = torch.cat([protein_feat, glycan_feat], dim=1)
        return self.regressor(combined)

# ==================== 训练流程 ====================
def train_model(model, train_loader, val_loader, epochs=100, lr=1e-4):
    model = model.to(device)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)
    
    best_val_loss = float('inf')
    train_losses, val_losses = [], []
    
    for epoch in range(epochs):
        model.train()
        train_loss = 0
        for protein, glycan, targets in train_loader:
            protein, glycan, targets = protein.to(device), glycan.to(device), targets.to(device)
            
            optimizer.zero_grad()
            outputs = model(protein, glycan)
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
        
        # 验证
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for protein, glycan, targets in val_loader:
                protein, glycan, targets = protein.to(device), glycan.to(device), targets.to(device)
                outputs = model(protein, glycan)
                val_loss += criterion(outputs, targets).item()
        
        # 记录损失
        train_loss /= len(train_loader)
        val_loss /= len(val_loader)
        train_losses.append(train_loss)
        val_losses.append(val_loss)
        
        # 保存最佳模型
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), "best_model.pth")
        
        print(f"Epoch {epoch+1}/{epochs} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")
    
    # 绘制损失曲线
    plt.plot(train_losses, label='Train Loss')
    plt.plot(val_losses, label='Val Loss')
    plt.xlabel('Epoch')
    plt.ylabel('MSE Loss')
    plt.legend()
    plt.show()
    
    return model

# ==================== 主流程 ====================
if __name__ == "__main__":
    # 模型训练配置
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    # 初始化模型
    model = TransformerRegressor(
        protein_dim=1280,
        glycan_dim=342,
        d_model=512,
        num_heads=8,
        num_layers=4,
        d_ff=2048
    )
    
    # 训练模型
    trained_model = train_model(model, train_loader, val_loader, epochs=10)
    
    # 测试评估
    trained_model.load_state_dict(torch.load("best_model.pth"))
    trained_model.eval()
    
    test_loss = 0
    criterion = nn.MSELoss()
    with torch.no_grad():
        for protein, glycan, targets in test_loader:
            protein, glycan, targets = protein.to(device), glycan.to(device), targets.to(device)
            outputs = model(protein, glycan)
            test_loss += criterion(outputs, targets).item()
    
    print(f"Test MSE Loss: {test_loss/len(test_loader):.4f}")