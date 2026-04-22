import csv

# 读取csv文件数据
data_Siglec = csv.reader(open("glycan_protein_716_Siglec.csv"))
data_lectin = csv.reader(open("glycan_protein_716_lectin.csv"))
data_ConA = csv.reader(open("glycan_protein_716_ConA.csv"))
data_Galectin = csv.reader(open("glycan_protein_716_Galectin.csv"))
data_train = csv.reader(open("glycan_protein_716_train.csv"))

# 定义氨基酸的单字母代码
amino_acids_codes = "ACDEFGHIKLMNPQRSTVWY"
# 生成氨基酸字典，使用标准的单字母代码
amino_acids_upper = {code: i + 1 for i, code in enumerate(amino_acids_codes)}
# 将大写字母字典转换为包含小写字母的字典
amino_acids = {**amino_acids_upper, **{k.lower(): v for k, v in amino_acids_upper.items()}}

# 最长的序列
max_sequence = 2710

time = 0

# 打开输出文件
with open("protein_encode_Siglec.csv", "w+", newline='') as f:
    writer = csv.writer(f)
    
    # 对氨基酸编号
    for row in data_Siglec:
        if time == 0:
            time += 1
            continue
        tmp = []
        for i in row[0]:
            if i.isdigit():
                continue
            tmp.append(amino_acids[i])
        zero = [0] * (max_sequence - len(tmp))
        tmp.extend(zero)
        # 去除第一列和最后一列的方括号
        writer.writerow(tmp)

time = 0

# 打开输出文件
with open("protein_encode_lectin.csv", "w+", newline='') as f:
    writer = csv.writer(f)
    
    # 对氨基酸编号
    for row in data_lectin:
        if time == 0:
            time += 1
            continue
        tmp = []
        for i in row[0]:
            if i.isdigit():
                continue
            tmp.append(amino_acids[i])
        zero = [0] * (max_sequence - len(tmp))
        tmp.extend(zero)
        # 去除第一列和最后一列的方括号
        writer.writerow(tmp)

time = 0

# 打开输出文件
with open("protein_encode_ConA.csv", "w+", newline='') as f:
    writer = csv.writer(f)
    
    # 对氨基酸编号
    for row in data_ConA:
        if time == 0:
            time += 1
            continue
        tmp = []
        for i in row[0]:
            if i.isdigit():
                continue
            tmp.append(amino_acids[i])
        zero = [0] * (max_sequence - len(tmp))
        tmp.extend(zero)
        # 去除第一列和最后一列的方括号
        writer.writerow(tmp)

time = 0

# 打开输出文件
with open("protein_encode_Galectin.csv", "w+", newline='') as f:
    writer = csv.writer(f)
    
    # 对氨基酸编号
    for row in data_Galectin:
        if time == 0:
            time += 1
            continue
        tmp = []
        for i in row[0]:
            if i.isdigit():
                continue
            tmp.append(amino_acids[i])
        zero = [0] * (max_sequence - len(tmp))
        tmp.extend(zero)
        # 去除第一列和最后一列的方括号
        writer.writerow(tmp)

time = 0

with open("protein_encode_train.csv", "w+", newline='') as f:
    writer = csv.writer(f)
    
    # 对氨基酸编号
    for row in data_train:
        if time == 0:
            time += 1
            continue
        tmp = []
        for i in row[0]:
            if i.isdigit():
                continue
            tmp.append(amino_acids[i])
        zero = [0] * (max_sequence - len(tmp))
        tmp.extend(zero)
        # 去除第一列和最后一列的方括号
        writer.writerow(tmp)