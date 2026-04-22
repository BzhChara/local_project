import csv
import itertools
import numpy as np

# 打开 CSV 文件并创建一个 reader 对象
with open("affinity_data.csv", newline='') as csvfile:
    reader = csv.reader(csvfile)
    
    # 初始化一个空列表来存储所有行
    matrix = []
    
    # 遍历 reader 对象，逐行读取数据
    for row in reader:
        # 将每行的每个字符串元素转换为适当的数据类型（例如 float），然后添加到矩阵列表
        matrix.append([value for value in row])

# 将空字符串替换为 None
for i in range(len(matrix)):
    matrix[i] = [None if value == '' else value for value in matrix[i]]

# 将其转换为 numpy 数组
matrix_array = np.array(matrix)

# print(matrix_array[0, 0])

result = []
shape = matrix_array.shape
f = open("数据整合.txt", "w+")

for i in range(1, shape[0]):
    for j in range(1, shape[1]):
        if matrix_array[i, j] is not None: 
            result.append(np.append(np.array(list(itertools.chain(*[matrix_array[i, 0].split(','), matrix_array[0, j].split(',')]))), matrix_array[i, j]))

for i in range(len(result)):
    f.write(str(result[i]).replace('\n','')+"\n")

f.close()
