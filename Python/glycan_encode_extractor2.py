import torch
import pandas as pd
import numpy as np
from torch import nn
import torch.nn.functional as F
from sklearn.preprocessing import StandardScaler

"""
# 固定随机种子
torch.manual_seed(42)
np.random.seed(42)
torch.cuda.manual_seed_all(42)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
"""

class GlycanAttentionExtractor(nn.Module):
    def __init__(self, input_dim=342, d_model=64):
        super().__init__()
        # 注意力机制参数
        self.W_query = nn.Linear(input_dim, d_model)
        self.W_key = nn.Linear(input_dim, d_model)
        self.W_value = nn.Linear(input_dim, d_model)
        
        # 特征融合层
        self.attention_fc = nn.Linear(input_dim + d_model, 32)
        
    def forward(self, x):
        # 注意力计算
        query = self.W_query(x)
        key = self.W_key(x)
        value = self.W_value(x)
        
        # 注意力权重
        scores = torch.matmul(query, key.transpose(1,2))
        attention_weights = F.normalize(scores, p=2, dim=1)
        
        # 加权特征
        weighted_values = torch.matmul(attention_weights, value)
        
        # 拼接原始特征
        glycan_feature = torch.cat((x, weighted_values), dim=1)
        
        # 最终特征投影
        return self.attention_fc(glycan_feature)

def process_glycan_features(input_file, output_file):
    # 数据加载与预处理
    df = pd.read_csv(input_file)
    glycan_data = df.values.astype(np.float32)
    
    # 标准化处理
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(glycan_data)
    
    # 转换为Tensor
    tensor_data = torch.FloatTensor(scaled_data)

    # 初始化模型
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = GlycanAttentionExtractor().to(device)

    # 特征提取
    model.eval()
    with torch.no_grad():
        features = model(tensor_data.to(device)).cpu().numpy()

    # 保存结果
    pd.DataFrame(features).to_csv(output_file, index=False, header=False)

if __name__ == "__main__":
    process_glycan_features(
        input_file="glycan_encode.csv",
        output_file="glycan_encode_new2.csv"
    )