import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import matplotlib.pyplot as plt

# 数据预处理
combined_data = pd.read_csv("combined_encode.csv", header=None)

protein_data = combined_data.iloc[:, :2710].values.astype(int)  
glycan_data = combined_data.iloc[:, 2710:2710+342].values.astype(float)  
targets = combined_data.iloc[:, -1].values.astype(float)

# 划分训练集和测试集（添加验证集）
X_prot_train, X_prot_test, X_gly_train, X_gly_test, y_train, y_test = train_test_split(
    protein_data, glycan_data, targets, test_size=0.2, random_state=42)
X_prot_train, X_prot_val, X_gly_train, X_gly_val, y_train, y_val = train_test_split(
    X_prot_train, X_gly_train, y_train, test_size=0.25, random_state=42)  

# 转换为PyTorch张量（修正糖链数据类型）
X_prot_train = torch.LongTensor(X_prot_train)
X_prot_val = torch.LongTensor(X_prot_val)
X_prot_test = torch.LongTensor(X_prot_test)
X_gly_train = torch.FloatTensor(X_gly_train)  
X_gly_val = torch.FloatTensor(X_gly_val)
X_gly_test = torch.FloatTensor(X_gly_test)
y_train = torch.FloatTensor(y_train)
y_val = torch.FloatTensor(y_val)
y_test = torch.FloatTensor(y_test)

# 数据加载器
batch_size = 32
train_dataset = torch.utils.data.TensorDataset(X_prot_train, X_gly_train, y_train)
val_dataset = torch.utils.data.TensorDataset(X_prot_val, X_gly_val, y_val)
test_dataset = torch.utils.data.TensorDataset(X_prot_test, X_gly_test, y_test)

train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=batch_size)
test_loader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size)

# 特征提取器
class ProteinFeatureExtractor(nn.Module):
    def __init__(self, vocab_size=21, embed_dim=64, hidden_dim=128, num_layers=2):
        super().__init__()
        self.embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=0)
        self.lstm = nn.LSTM(embed_dim, hidden_dim, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_dim, 64)

    def forward(self, x):
        x = self.embedding(x)  # (batch, seq_len, embed_dim)
        _, (h_n, _) = self.lstm(x)  # h_n: (num_layers, batch, hidden_dim)
        return self.fc(h_n[-1])  # (batch, 64)

class GlycanFeatureExtractor(nn.Module):
    def __init__(self, input_dim=342, hidden_dim=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.ReLU(),
            nn.Linear(256, hidden_dim)
        )

    def forward(self, x):
        return self.net(x)  # (batch, hidden_dim)

# Transformer模型
class RegressionTransformer(nn.Module):
    def __init__(self, d_model=64, nhead=8, num_layers=4):
        super().__init__()
        # 蛋白质特征提取器添加残差连接
        self.prot_extractor = nn.Sequential(
            ProteinFeatureExtractor(),
            nn.LayerNorm(64)
        )
        # 糖链特征提取器加深
        self.gly_extractor = nn.Sequential(
            GlycanFeatureExtractor(),
            nn.LayerNorm(64),
            nn.Linear(64, 128),
            nn.ReLU(),
            nn.Linear(128, 64)
        )
        
        # Transformer添加跳跃连接
        self.transformer = nn.TransformerEncoder(
            encoder_layer=nn.TransformerEncoderLayer(
                d_model=d_model,
                nhead=nhead,
                dim_feedforward=256,
                batch_first=True,
                dropout=0.1
            ),
            num_layers=num_layers
        )
        self.transformer_norm = nn.LayerNorm(d_model)

        self.regressor = nn.Sequential(
            nn.Linear(d_model, 32),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(32, 1)
        )
        
        self.pos_encoder = nn.Parameter(torch.randn(1, 2, d_model))  # 可学习位置编码

    def forward(self, prot, gly):
        # 特征提取
        prot_feat = self.prot_extractor(prot)  # (batch, 64)
        gly_feat = self.gly_extractor(gly)     # (batch, 64)
        
        # 合并特征
        combined = torch.stack([prot_feat, gly_feat], dim=1)  # (batch, 2, 64)
        combined += self.pos_encoder
        
        # Transformer处理
        transformer_out = self.transformer(combined)
        transformer_out = self.transformer_norm(transformer_out) # 添加归一化
        
        # 回归预测
        return self.regressor(transformer_out.mean(dim=1)).squeeze()

# 模型初始化
model = RegressionTransformer()

# 权重初始化
def init_weights(m):
    if isinstance(m, nn.Linear):
        nn.init.kaiming_normal_(m.weight, nonlinearity='relu')
        if m.bias is not None:
            nn.init.constant_(m.bias, 0)
    elif isinstance(m, nn.Embedding):
        nn.init.normal_(m.weight, mean=0, std=0.1)
model.apply(init_weights)

# 训练配置
criterion = nn.MSELoss()
optimizer = optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=5, factor=0.5)
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model.to(device)

# 训练循环
best_val_loss = float('inf')
patience = 10
epochs = 10

for epoch in range(epochs):
    # 训练阶段
    model.train()
    train_loss = 0
    for prot, gly, y in train_loader:
        prot, gly, y = prot.to(device), gly.to(device), y.to(device)
        
        optimizer.zero_grad()
        outputs = model(prot, gly)
        loss = criterion(outputs, y)
        loss.backward()

        """""
        # 检查梯度是否存在
        total_grad = 0
        for name, param in model.named_parameters():
            if param.grad is not None:
                total_grad += param.grad.abs().sum().item()
        print(f"Batch gradient magnitude: {total_grad:.4f}")
        """""

        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        
        train_loss += loss.item() * prot.size(0)

    # 验证阶段
    model.eval()
    val_loss = 0
    with torch.no_grad():
        for prot, gly, y in val_loader:
            prot, gly, y = prot.to(device), gly.to(device), y.to(device)
            outputs = model(prot, gly)
            val_loss += criterion(outputs, y).item() * prot.size(0)
    
    # 计算平均损失
    train_loss /= len(train_loader.dataset)
    val_loss /= len(val_loader.dataset)
    scheduler.step(val_loss)
    
    # 早停检查
    if val_loss < best_val_loss:
        best_val_loss = val_loss
        torch.save(model.state_dict(), 'best_model.pth')
        patience_counter = 0
    else:
        patience_counter += 1
        if patience_counter >= patience:
            print(f'Early stopping at epoch {epoch+1}')
            break
    
    print(f'Epoch {epoch+1:3d} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}')

# 加载最佳模型
model.load_state_dict(torch.load('best_model.pth'))

# 测试评估
model.eval()
y_true, y_pred = [], []
with torch.no_grad():
    for prot, gly, y in test_loader:
        prot, gly = prot.to(device), gly.to(device)
        outputs = model(prot, gly).cpu().numpy()
        y_true.extend(y.cpu().numpy())
        y_pred.extend(outputs.cpu().numpy())

y_true = np.array(y_true)
y_pred = np.array(y_pred)

# 计算指标
mae = mean_absolute_error(y_true, y_pred)
mse = mean_squared_error(y_true, y_pred)
r2 = r2_score(y_true, y_pred)
print(f'MAE: {mae:.4f}')
print(f'MSE: {mse:.4f}')
print(f'R²:  {r2:.4f}')

# 残差图
plt.figure(figsize=(10,6))
plt.scatter(y_true, y_pred - y_true, alpha=0.5, edgecolors='w')
plt.hlines(0, y_true.min(), y_true.max(), colors='red', linestyles='dashed')
plt.xlabel('True Values', fontsize=12)
plt.ylabel('Residuals', fontsize=12)
plt.title('Regression Residual Plot', fontsize=14)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()