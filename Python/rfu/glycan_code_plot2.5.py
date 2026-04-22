import re
import matplotlib.pyplot as plt
from collections import Counter

def extract_main_chain_lengths(filename):
    # 提取主链长度
    chain_lengths = []
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
            main_chain = re.sub(r'\([^)]*\)', '', glycan_chain)
            sugar_names = re.findall(r'[A-Za-z]+', main_chain)
            chain_length = len(sugar_names)
            chain_lengths.append(chain_length)
    return chain_lengths

def plot_main_chain_length_distribution(chain_lengths):
    # 统计每种主链长度的频次
    length_counts = Counter(chain_lengths)
    
    # 提取数据用于绘图
    lengths = list(length_counts.keys())
    counts = list(length_counts.values())
    
    plt.figure(figsize=(16, 8))
    bars = plt.bar(lengths, counts, color='#2c7bb6', edgecolor='white', alpha=0.8, width=0.5)
    
    # 设置刻度与标签样式
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('Main Chain Length', fontsize=14)
    plt.ylabel('Number of Oligosaccharides', fontsize=14)
    
    # 确保x轴刻度完全显示
    plt.xticks(lengths, fontsize=12)
    
    # 添加数据标签
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{height}',
                 ha='center', va='bottom', fontsize=12)
    
    # 添加网格线和标题
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.title('Main Chain Length Distribution', fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    plt.show()

# 示例用法
if __name__ == "__main__":
    input_file = 'sugarlist.txt'  # 替换为你的输入文件路径
    chain_lengths = extract_main_chain_lengths(input_file)
    plot_main_chain_length_distribution(chain_lengths)