import pandas as pd
from scipy.stats import zscore

# 读取CSV文件 第一行默认为表头
df = pd.read_csv('glycan_protein_716.csv')

# 提取第一列
first_column = df.iloc[:, 0]  # 第一列数据

# 对每一行的数据进行 Z-score 标准化
data_to_standardize = df.iloc[:, 2:]
zscore_data = data_to_standardize.apply(lambda x: zscore(x, nan_policy='omit'), axis=1)

# 添加第一列数据到标准化数据
zscore_data.insert(0, first_column.name, first_column)

# 将处理后的数据输出到新的CSV文件
zscore_data.to_csv('z-output.csv', index=False)




