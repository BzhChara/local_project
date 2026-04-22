import re
import csv

# 读取csv文件数据
data = csv.reader(open("glycan_protein_716.csv", encoding="utf-8"))

# 训练集
train = []

time = 0

a = open("glycan_protein_716_Siglec.csv", "w+", newline='')
b = open("glycan_protein_716_lectin.csv", "w+", newline='')
c = open("glycan_protein_716_ConA.csv", "w+", newline='')
d = open("glycan_protein_716_Galectin.csv", "w+", newline='')

writer_a = csv.writer(a)
writer_b = csv.writer(b)
writer_c = csv.writer(c)
writer_d = csv.writer(d)


for row in data:
    if time == 0:
        time += 1
        head = row[2:]
        writer_a.writerow(head)
        writer_b.writerow(head)
        writer_c.writerow(head)
        writer_d.writerow(head)
        continue
    if re.match("Siglec", row[1]):
        writer_a.writerow(row[2:])
    elif re.match("C-type lectin", row[1]):
        writer_b.writerow(row[2:])
    elif re.match("ConA", row[1]):
        writer_c.writerow(row[2:])
    elif re.match("Galectin", row[1]) and re.match("Human", row[0]):
        writer_d.writerow(row[2:])
    else:
        train.append(row[2:])

a.close()
b.close()
c.close()
d.close()

# 输出训练集
with open("glycan_protein_716_train.csv", "w+", newline='') as f:
    writer = csv.writer(f)

    writer.writerow(head)
    for i in train:
        writer.writerow(i)