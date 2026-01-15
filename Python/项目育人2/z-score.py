import pandas as pd
from scipy.stats import zscore

# 读取CSV文件 第一行默认为表头
df = pd.read_csv('glycan_protein_716_Siglec.csv')

# 提取第一列
third_column = df.iloc[:, 0] 

# 对每一行的数据进行 Z-score 标准化
data_to_standardize = df.iloc[:, 1:]
zscore_data = data_to_standardize.apply(lambda x: zscore(x, nan_policy='omit'), axis=1)

# 添加第1,2,3列数据到标准化数据
zscore_data.insert(0, third_column.name, third_column)

# 将处理后的数据输出到新的CSV文件
zscore_data.to_csv('z-output_Siglec.csv', index=False)

# 读取CSV文件 第一行默认为表头
df = pd.read_csv('glycan_protein_716_lectin.csv')

# 提取第一列
third_column = df.iloc[:, 0] 

# 对每一行的数据进行 Z-score 标准化
data_to_standardize = df.iloc[:, 1:]
zscore_data = data_to_standardize.apply(lambda x: zscore(x, nan_policy='omit'), axis=1)

# 添加第1,2,3列数据到标准化数据
zscore_data.insert(0, third_column.name, third_column)

# 将处理后的数据输出到新的CSV文件
zscore_data.to_csv('z-output_lectin.csv', index=False)

# 读取CSV文件 第一行默认为表头
df = pd.read_csv('glycan_protein_716_ConA.csv')

# 提取第一列
third_column = df.iloc[:, 0] 

# 对每一行的数据进行 Z-score 标准化
data_to_standardize = df.iloc[:, 1:]
zscore_data = data_to_standardize.apply(lambda x: zscore(x, nan_policy='omit'), axis=1)

# 添加第1,2,3列数据到标准化数据
zscore_data.insert(0, third_column.name, third_column)

# 将处理后的数据输出到新的CSV文件
zscore_data.to_csv('z-output_ConA.csv', index=False)

# 读取CSV文件 第一行默认为表头
df = pd.read_csv('glycan_protein_716_Galectin.csv')

# 提取第一列
third_column = df.iloc[:, 0] 

# 对每一行的数据进行 Z-score 标准化
data_to_standardize = df.iloc[:, 1:]
zscore_data = data_to_standardize.apply(lambda x: zscore(x, nan_policy='omit'), axis=1)

# 添加第1,2,3列数据到标准化数据
zscore_data.insert(0, third_column.name, third_column)

# 将处理后的数据输出到新的CSV文件
zscore_data.to_csv('z-output_Galectin.csv', index=False)

# 读取CSV文件 第一行默认为表头
df = pd.read_csv('glycan_protein_716_train.csv')

# 提取第一列
third_column = df.iloc[:, 0] 

# 对每一行的数据进行 Z-score 标准化
data_to_standardize = df.iloc[:, 1:]
zscore_data = data_to_standardize.apply(lambda x: zscore(x, nan_policy='omit'), axis=1)

# 添加第1,2,3列数据到标准化数据
zscore_data.insert(0, third_column.name, third_column)

# 将处理后的数据输出到新的CSV文件
zscore_data.to_csv('z-output_train.csv', index=False)