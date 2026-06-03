import re
import matplotlib.pyplot as plt
from collections import Counter

# 读取糖链数据并统计支链数量
def count_branches(glycan_chain):
    branch_count = len(re.findall(r'\([^)]*\)', glycan_chain))
    return branch_count

def extract_branch_counts(filename):
    branch_counts = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 3):
            if i + 2 >= len(lines):
                break
            single_line = lines[i].strip()
            if not single_line.startswith('single;'):
                continue
            parts = single_line.split(';')
            if len(parts) < 2:
                continue
            glycan_chain = parts[1]
            branch_count = count_branches(glycan_chain)
            branch_counts.append(branch_count)
    return branch_counts

# 绘制柱状图
def plot_branch_count_distribution(branch_counts):
    branch_count_dist = Counter(branch_counts)
    branches = list(branch_count_dist.keys())
    counts = list(branch_count_dist.values())

    plt.figure(figsize=(16, 8))
    bars = plt.bar(branches, counts, color='#2c7bb6', edgecolor='white', alpha=0.8, width=0.5)

    # 设置刻度与标签样式
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('Number of Branches', fontsize=14)
    plt.ylabel('Number of Oligosaccharides', fontsize=14)

    # 添加数据标签
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{height}',
                 ha='center', va='bottom', fontsize=12)

    # 添加网格线和标题
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.title('Glycan Chain Branch Count Distribution', fontsize=16, fontweight='bold')

    plt.tight_layout()
    plt.show()

# 示例用法
if __name__ == "__main__":
    import numpy as np
    input_file = 'sugarlist.txt'  # 替换为你的输入文件路径
    branch_counts = extract_branch_counts(input_file)
    plot_branch_count_distribution(branch_counts)