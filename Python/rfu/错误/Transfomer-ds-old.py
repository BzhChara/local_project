import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import matplotlib.pyplot as plt

# 读取数据
combined_data = pd.read_csv("combined_encode.csv", header=None)
protein_data = combined_data.iloc[:, :2710].values.astype(int)
glycan_data = combined_data.iloc[:, 2710:2710+342].values.astype(int)
targets = combined_data.iloc[:, -1].values.astype(float)

# 划分训练集和测试集
X_prot_train, X_prot_test, X_gly_train, X_gly_test, y_train, y_test = train_test_split(
    protein_data, glycan_data, targets, test_size=0.2, random_state=42
)

# 转换为PyTorch张量
X_prot_train = torch.LongTensor(X_prot_train)
X_prot_test = torch.LongTensor(X_prot_test)
X_gly_train = torch.LongTensor(X_gly_train)
X_gly_test = torch.LongTensor(X_gly_test)
y_train = torch.FloatTensor(y_train)
y_test = torch.FloatTensor(y_test)

# 定义数据加载器
batch_size = 32
train_dataset = torch.utils.data.TensorDataset(X_prot_train, X_gly_train, y_train)
train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, shuffle=True)

# 定义特征提取模型
class ProteinFeatureExtractor(nn.Module):
    def __init__(self, vocab_size=21, embedding_dim=64, hidden_dim=128, num_layers=2):
        super().__init__()
        self.embedding = nn.Embedding(vocab_size, embedding_dim)
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_dim, 64)

    def forward(self, x):
        x = self.embedding(x)
        lstm_out, _ = self.lstm(x)
        return self.fc(lstm_out[:, -1, :])

class GlycanFeatureExtractor(nn.Module):
    def __init__(self, input_dim=342, embedding_dim=64):
        super().__init__()
        self.fc = nn.Linear(input_dim, embedding_dim)
        self.query = nn.Linear(embedding_dim, embedding_dim)
        self.key = nn.Linear(embedding_dim, embedding_dim)
        self.value = nn.Linear(embedding_dim, embedding_dim)
        self.fc_out = nn.Linear(embedding_dim, 64)

    def forward(self, x):
        x = self.fc(x.float())
        x = x.unsqueeze(1)
        Q = self.query(x)
        K = self.key(x)
        V = self.value(x)
        attn = F.softmax(torch.matmul(Q, K.transpose(1,2)) / np.sqrt(64), dim=-1)
        out = torch.matmul(attn, V).squeeze(1)
        return self.fc_out(out)

class PositionalEncoding(nn.Module):
    def __init__(self, d_model=64, max_len=2):
        super().__init__()
        self.position_embedding = nn.Parameter(torch.randn(1, max_len, d_model))
    
    def forward(self, x):
        # x shape: [batch_size, seq_len, d_model]
        return x + self.position_embedding[:, :x.size(1), :]

class TransformerRegressionModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.protein_extractor = ProteinFeatureExtractor()
        self.glycan_extractor = GlycanFeatureExtractor()
        self.transformer = nn.TransformerEncoder(
            encoder_layer=nn.TransformerEncoderLayer(
                d_model=64,  # 必须与特征维度匹配
                nhead=8,
                dim_feedforward=256,
                batch_first=True  # 添加batch_first参数
            ),
            num_layers=3
        )
        self.regressor = nn.Sequential(
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 1)
        )
        self.pos_encoder = PositionalEncoding(d_model=64)

    def forward(self, prot, gly):
        prot_feat = self.protein_extractor(prot)
        gly_feat = self.glycan_extractor(gly)
        combined = torch.stack([prot_feat, gly_feat], dim=1)
        combined = self.pos_encoder(combined)
        transformer_out = self.transformer(combined)
        return self.regressor(transformer_out.mean(dim=1)).squeeze()

# 训练配置
model = TransformerRegressionModel()
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)
epochs = 5

# 训练循环
for epoch in range(epochs):
    model.train()
    total_loss = 0
    for prot, gly, y in train_loader:
        optimizer.zero_grad()
        outputs = model(prot, gly)
        loss = criterion(outputs, y)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    print(f"Epoch {epoch+1}, Loss: {total_loss/len(train_loader):.4f}")

# 评估模型
model.eval()
with torch.no_grad():
    y_pred = model(X_prot_test, X_gly_test)
    y_pred = y_pred.numpy()
    y_true = y_test.numpy()

# 计算指标
mae = mean_absolute_error(y_true, y_pred)
mse = mean_squared_error(y_true, y_pred)
r2 = r2_score(y_true, y_pred)

print(f"MAE: {mae:.4f}")
print(f"MSE: {mse:.4f}")
print(f"R²: {r2:.4f}")

# 绘制残差图
plt.figure(figsize=(8,6))
plt.scatter(y_true, y_pred - y_true, alpha=0.5)
plt.hlines(0, y_true.min(), y_true.max(), colors='red')
plt.xlabel("True Values")
plt.ylabel("Residuals")
plt.title("Residual Plot")
plt.show()