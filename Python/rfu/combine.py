import csv
import pandas as pd
# 读取蛋白质文件
protein_data = pd.read_csv("protein_esm.csv", header=None)
# 读取寡糖编码文件
glycan_data = pd.read_csv("glycan_encode.csv", header=None)
# 读取值文件
values_data = pd.read_csv("glycan_protein_rfu.csv", header=None)

# 打开输出文件
with open("combined_encode_esm.csv", "w+", newline='') as f:
    writer = csv.writer(f)
    
    # 遍历蛋白质
    for i, protein_row in protein_data.iterrows():
        # 遍历寡糖编码
        for j, glycan_row in glycan_data.iterrows():
            # 检查值文件中的对应值是否为空白
            value = values_data.iloc[i + 1, j + 1]
            if pd.isna(value):
                continue
            # 合并蛋白质和寡糖编码
            combined_row = list(protein_row) + list(glycan_row)
            # 加上值文件中的值
            combined_row.append(float(value))
            # 写入输出文件
            writer.writerow(combined_row)

