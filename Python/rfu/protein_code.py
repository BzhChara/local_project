import csv

# 读取csv文件数据
data = csv.reader(open("glycan_protein_rfu.csv"))

# 定义氨基酸的单字母代码
amino_acids_codes = "ACDEFGHIKLMNPQRSTVWY"
# 生成氨基酸字典，使用标准的单字母代码
amino_acids_upper = {code: i + 1 for i, code in enumerate(amino_acids_codes)}
# 将大写字母字典转换为包含小写字母的字典
amino_acids = {**amino_acids_upper, **{k.lower(): v for k, v in amino_acids_upper.items()}}

time = 0
num = []

# 统计最长的序列
for row in data:
    if time == 0:
            time += 1
            continue
    num.append(sum(1 for char in row[0] if char.isalpha()))
max_sequence = max(num)
min_sequence = min(num)
print(max_sequence, min_sequence)

time = 0

# 重新读取csv文件数据
data = csv.reader(open("glycan_protein_rfu.csv"))

# 打开输出文件
with open("protein_encode.csv", "w+", newline='') as f:
    writer = csv.writer(f)
    
    # 对氨基酸编号
    for row in data:
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
