from __future__ import print_function
import csv
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

input_size = 2710  # 每个特征数字列表的长度
hidden_size = 32
num_layers = 2

lstm = nn.LSTM(input_size, hidden_size, num_layers, batch_first=True)

sequence_data = []
with open("蛋白质编码.txt", "r") as file:
    for line in file:
        line = line[1: -2]
        features = list(map(float, line.strip().split(',')))
        sequence_tensor = torch.tensor(features, dtype=torch.float32)
        sequence_data.append(sequence_tensor)

# 将序列数据转换为批次数据
sequence_tensors = torch.stack(sequence_data, dim=0) 

# 定义批次大小
batch_size = 32

# 计算批次数量
num_batches = (len(sequence_tensors) + batch_size - 1) // batch_size

# 处理数据并获取LSTM层的输出
all_outputs = []  # 存储所有批次的输出
for batch_idx in range(num_batches):
    # 计算批次的起始和结束索引
    start_idx = batch_idx * batch_size
    end_idx = start_idx + batch_size
    # 获取当前批次的数据
    batch_data = sequence_tensors[start_idx:end_idx]
    
    # 将批次数据传递给LSTM层，PyTorch将自动初始化隐藏状态和细胞状态
    outputs, _ = lstm(batch_data)  # 注意这里我们不需要提供隐藏和细胞状态作为输入

    # 存储当前批次的输出
    all_outputs.append(outputs)

# 将所有批次的输出合并
final_output = torch.cat(all_outputs, dim=0)

# 将Tensor转换为NumPy数组
final_output_array = final_output.detach().cpu().numpy()

# 将NumPy数组保存到TXT文件中
with open("蛋白质特征提取.txt", "w") as f:
    for row in final_output_array:
        # 将数组的每个元素转换为字符串，并使用逗号分隔
        row_str = ','.join(map(str, row))
        # 写入字符串到文件，并添加换行符
        f.write("[" + row_str + "]" + '\n')




        