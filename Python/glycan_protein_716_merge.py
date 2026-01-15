import pandas as pd

# 读取三个CSV文件
original_df = pd.read_csv('glycan_protein_716.csv')
protein_encode_df = pd.read_csv('protein_esm.csv', header=None)
glycan_encode_df = pd.read_csv('glycan_encode.csv', header=None)

# 将蛋白质编码替换到原始文件中
original_df['protein sequence'] = protein_encode_df.apply(lambda row: row.values.tolist(), axis=1)

# 保存为新的CSV文件
original_df.to_csv('glycan_protein_716_new.csv', index=False)

