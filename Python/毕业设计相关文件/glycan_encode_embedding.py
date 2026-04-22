import torch
import pandas as pd
import numpy as np
from torch import nn
import torch.nn.functional as F
from sklearn.preprocessing import StandardScaler

glycan = pd.read_csv("glycan_encode.csv",header=None)
max_value = glycan.max().max()
glycan_data = glycan.values.tolist()
embedding = nn.Embedding(max_value +1,1)
input = torch.LongTensor(glycan_data)
e = embedding(input)
feature = []
for i in e:
    j = i.view(1,-1)
    feature.append(j)
my_feature =[]
for m in feature:
    n = m.tolist()
    for k in n:
        k = list(k)
        my_feature.append(k)
my_feature = np.array(my_feature)
glycan_feature = pd.DataFrame(my_feature)

# 保存结果
pd.DataFrame(glycan_feature).to_csv("glycan_encode_embedding.csv", index=False, header=False)

