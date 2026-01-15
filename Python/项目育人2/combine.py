import csv
import pandas as pd
# 读取蛋白质编码文件
protein_data_Siglec = pd.read_csv("protein_encode_Siglec.csv", header=None)
protein_data_lectin = pd.read_csv("protein_encode_lectin.csv", header=None)
protein_data_ConA = pd.read_csv("protein_encode_ConA.csv", header=None)
protein_data_Galectin = pd.read_csv("protein_encode_Galectin.csv", header=None)
protein_data_train = pd.read_csv("protein_encode_train.csv", header=None)
# 读取寡糖编码文件
glycan_data = pd.read_csv("glycan_encode.csv", header=None)
# 读取值文件
values_data_Siglec = pd.read_csv("z-output_Siglec.csv", header=None)
values_data_lectin = pd.read_csv("z-output_lectin.csv", header=None)
values_data_ConA = pd.read_csv("z-output_ConA.csv", header=None)
values_data_Galectin = pd.read_csv("z-output_Galectin.csv", header=None)
values_data_train = pd.read_csv("z-output_train.csv", header=None)

# 打开输出文件
with open("combined_encode_Siglec.csv", "w+", newline='') as f:
    writer = csv.writer(f)
    
    # 遍历蛋白质编码
    for i, protein_row in protein_data_Siglec.iterrows():
        # 遍历寡糖编码
        for j, glycan_row in glycan_data.iterrows():
            # 检查值文件中的对应值是否为空白
            value = values_data_Siglec.iloc[i + 1, j + 1]
            if pd.isna(value):
                continue
            # 合并蛋白质和寡糖编码
            combined_row = list(protein_row) + list(glycan_row)
            # 加上值文件中的值
            combined_row.append(float(value))
            # 写入输出文件
            writer.writerow(combined_row)

# 打开输出文件
with open("combined_encode_lectin.csv", "w+", newline='') as f:
    writer = csv.writer(f)
    
    # 遍历蛋白质编码
    for i, protein_row in protein_data_lectin.iterrows():
        # 遍历寡糖编码
        for j, glycan_row in glycan_data.iterrows():
            # 检查值文件中的对应值是否为空白
            value = values_data_lectin.iloc[i + 1, j + 1]
            if pd.isna(value):
                continue
            # 合并蛋白质和寡糖编码
            combined_row = list(protein_row) + list(glycan_row)
            # 加上值文件中的值
            combined_row.append(float(value))
            # 写入输出文件
            writer.writerow(combined_row)

# 打开输出文件
with open("combined_encode_ConA.csv", "w+", newline='') as f:
    writer = csv.writer(f)
    
    # 遍历蛋白质编码
    for i, protein_row in protein_data_ConA.iterrows():
        # 遍历寡糖编码
        for j, glycan_row in glycan_data.iterrows():
            # 检查值文件中的对应值是否为空白
            value = values_data_ConA.iloc[i + 1, j + 1]
            if pd.isna(value):
                continue
            # 合并蛋白质和寡糖编码
            combined_row = list(protein_row) + list(glycan_row)
            # 加上值文件中的值
            combined_row.append(float(value))
            # 写入输出文件
            writer.writerow(combined_row)

            # 打开输出文件
with open("combined_encode_Galectin.csv", "w+", newline='') as f:
    writer = csv.writer(f)
    
    # 遍历蛋白质编码
    for i, protein_row in protein_data_Galectin.iterrows():
        # 遍历寡糖编码
        for j, glycan_row in glycan_data.iterrows():
            # 检查值文件中的对应值是否为空白
            value = values_data_Galectin.iloc[i + 1, j + 1]
            if pd.isna(value):
                continue
            # 合并蛋白质和寡糖编码
            combined_row = list(protein_row) + list(glycan_row)
            # 加上值文件中的值
            combined_row.append(float(value))
            # 写入输出文件
            writer.writerow(combined_row)

# 打开输出文件
with open("combined_encode_train.csv", "w+", newline='') as f:
    writer = csv.writer(f)
    
    # 遍历蛋白质编码
    for i, protein_row in protein_data_train.iterrows():
        # 遍历寡糖编码
        for j, glycan_row in glycan_data.iterrows():
            # 检查值文件中的对应值是否为空白
            value = values_data_train.iloc[i + 1, j + 1]
            if pd.isna(value):
                continue
            # 合并蛋白质和寡糖编码
            combined_row = list(protein_row) + list(glycan_row)
            # 加上值文件中的值
            combined_row.append(float(value))
            # 写入输出文件
            writer.writerow(combined_row)