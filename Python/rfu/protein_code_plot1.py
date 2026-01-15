import csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 读取蛋白质序列数据
with open('glycan_protein_rfu.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)  # 跳过标题行
    proteins = [row[0] for row in reader]

# 计算蛋白质长度
protein_lengths = [len(seq) for seq in proteins]

# --------------------------
# 直方图 (Histogram)
# --------------------------
plt.figure(figsize=(16, 8))
max_len = max(protein_lengths)
end = ((max_len // 50) + 1) * 50  
bins = np.arange(0, end + 1, 50)  

# 绘制直方图（直接使用bins边界）
plt.hist(protein_lengths, 
         bins=bins,          # 强制对齐刻度
         color='#2c7bb6', 
         edgecolor='white',
         alpha=0.8)

# --------------------------
# 设置刻度与标签样式
# --------------------------
xticks = np.arange(0, end + 1, 50)  
plt.xticks(xticks, 
           rotation=45,       # 刻度倾斜45度
           ha='right',        # 标签右对齐（与旋转配合防重叠）
           fontsize=12)
plt.xlim(0, end)             # 固定坐标轴范围

plt.title('Protein Sequence Length Distribution Histogram', fontsize=16, fontweight='bold')
plt.xlabel('Sequence Length', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.grid(axis='y', linestyle='--', alpha=0.7)

# 添加数据统计标注
stats_text = f"""Total Proteins: {len(protein_lengths)}
Mean Length: {np.mean(protein_lengths):.1f}
Max Length: {max(protein_lengths)}"""
plt.annotate(stats_text, 
             xy=(0.95, 0.95), 
             xycoords='axes fraction',
             ha='right', 
             va='top',
             bbox=dict(boxstyle='round', alpha=0.8, facecolor='white'))

plt.tight_layout()
plt.show()

# --------------------------
# 小提琴图 (Violin Plot)
# --------------------------
plt.figure(figsize=(16, 8))
sns.violinplot(y=protein_lengths,
               inner='quartile',  # 显示四分位线
               palette='Blues',
               linewidth=1)

plt.title('Protein Length Distribution Violin Plot', fontsize=16, fontweight='bold')
plt.ylabel('Sequence Length', fontsize=14)
plt.grid(axis='y', linestyle='--', alpha=0.7)

# 添加数值标记
q1, median, q3 = np.percentile(protein_lengths, [25, 50, 75])
plt.text(0.3, median+5, f'Median: {median}', 
         color='darkred', 
         fontweight='bold')
plt.text(0.3, q1-5, f'Q1: {q1}', 
         color='navy', 
         verticalalignment='top')
plt.text(0.3, q3+5, f'Q3: {q3}', 
         color='navy', 
         verticalalignment='bottom')

plt.tight_layout()
plt.show()