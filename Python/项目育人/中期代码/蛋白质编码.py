import csv

#读取csv文件数据
data=csv.reader(open("glycan_protein_rfu.csv"))
#print(data)

#建立字典
# 定义氨基酸的单字母代码
amino_acids_codes = "ACDEFGHIKLMNPQRSTVWY"

# 生成氨基酸字典，使用标准的单字母代码
amino_acids_upper = {code: i + 1 for i, code in enumerate(amino_acids_codes)}


# 将大写字母字典转换为包含小写字母的字典
amino_acids = {**amino_acids_upper, **{k.lower(): v for k, v in amino_acids_upper.items()}}

time=0
num=[]
f=open("蛋白质分析.txt","w+")

#统计最长的序列
for row in data:
    num.append(sum(1 for char in row[0] if char.isalpha()))
max_sequence=max(num)
print(max_sequence)

#读取csv文件数据
data=csv.reader(open("glycan_protein_rfu.csv"))
#print(data)

#对氨基酸编号
for row in data:
    if time==0:
        time=time+1
        continue
    tmp=[]
    for i in row[0]:
        if i.isdigit():
            continue
        tmp.append(amino_acids[i])
        #print(amino_acids[i])
    zero=[0]*(max_sequence-len(tmp))
    tmp.extend(zero)
    f.write(str(tmp)+"\n")

f.close()