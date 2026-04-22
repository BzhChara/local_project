import csv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import tab20

# 标准氨基酸列表 (按字母顺序)
amino_acids = sorted(['A','C','D','E','F','G','H','I','K','L',
                     'M','N','P','Q','R','S','T','V','W','Y'])

# 读取数据
with open('glycan_protein_rfu.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)  # 跳过标题
    proteins = [row[0] for row in reader]

# 初始化存储矩阵
protein_data = []

# 统计每个蛋白质的氨基酸组成
for seq in proteins:
    counts = {aa:0 for aa in amino_acids}
    total = 0
    for char in seq:
        if char in counts:
            counts[char] += 1
            total += 1
    
    percentages = []
    for aa in amino_acids:
        percentages.append(counts[aa]/total*100 if total !=0 else 0)
    
    protein_data.append(percentages)

# 创建绘图数据
plt.figure(figsize=(16, 8))
ax = plt.gca()
x = np.arange(len(proteins))  # 蛋白质编号
colors = tab20.colors  # 颜色方案

# 创建堆叠柱状图
bottom = np.zeros(len(proteins))
for idx, aa in enumerate(amino_acids):
    values = [protein_data[i][idx] for i in range(len(proteins))]
    ax.bar(x, values, bottom=bottom, label=aa, color=colors[idx])
    bottom += values

# 设置坐标轴格式
max_x = len(proteins) - 1
xticks_step = 10
xticks_pos = np.arange(0, max_x + xticks_step, xticks_step)
ax.set_xticks(xticks_pos)
ax.set_xticklabels(xticks_pos, rotation=45, ha='right')

ax.set_title('Amino Acid Composition per Protein', fontsize=16, fontweight='bold')
ax.set_xlabel('Protein Index', fontsize=14)
ax.set_ylabel('Percentage (%)', fontsize=14)
ax.set_ylim(0, 100)
ax.grid(axis='y', linestyle='--', alpha=0.7)

# 图例优化
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], 
         bbox_to_anchor=(1.12, 1),
         loc='upper left',
         ncol=2,
         title='Amino Acids',
         title_fontsize=10)

plt.tight_layout()
plt.show()