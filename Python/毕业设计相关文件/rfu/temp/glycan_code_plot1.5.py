import re
import matplotlib.pyplot as plt
from collections import Counter

def read_triplets(file_path):
    """生成器函数，每次读取三行作为一组"""
    with open(file_path, 'r') as file:
        lines = []
        for line in file:
            lines.append(line.strip())
            if len(lines) == 3:
                yield lines
                lines = []

def extract_sugar_counts(line):
    """精确提取糖类型数量"""
    counts_part = line.split(';')[-1]
    numbers = re.findall(r':(\d+)', counts_part)
    return sum(map(int, numbers))

def plot_sugar_distribution(glycan_data):
    """绘制水平标签柱状图"""
    # 数据收集
    single = []
    double = []
    triple = []
    
    for triplet in glycan_data:
        single.append(extract_sugar_counts(triplet[0]))
        double.append(extract_sugar_counts(triplet[1]))
        triple.append(extract_sugar_counts(triplet[2]))
    
    # 统计频率
    single_dist = Counter(single)
    double_dist = Counter(double)
    triple_dist = Counter(triple)

    # 绘图设置
    plt.figure(figsize=(18, 8))  # 适当加宽画布
    colors = ['blue', 'red', 'green']
    
    # 创建子图
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    plt.subplots_adjust(wspace=0.35)  # 增加子图间距

    # 通用绘图函数
    def create_plot(ax, data, color, title):
        keys = sorted(data.keys())
        values = [data[k] for k in keys]
        
        # 柱状图参数
        ax.bar(keys, values, 
               color=color, 
               alpha=0.7, 
               width=0.8,
               edgecolor='black',
               linewidth=0.5)
        
        # 智能标签间隔
        max_visible = 15  # 最大可见标签数
        step = max(1, len(keys) // max_visible)
        visible_ticks = keys[::step]
        
        # 设置坐标轴
        ax.set_xticks(visible_ticks)
        ax.set_xticklabels(visible_ticks, 
                          fontsize=10,  # 保持原字号
                          rotation=0,    # 关键修改：去除旋转
                          ha='center')   # 居中对齐
        
        # 设置标签和标题
        ax.set_title(title, 
                    fontsize=16,
                    fontweight='bold',
                    pad=12)
        ax.set_xlabel('Sugar Units', 
                     fontsize=14)
        ax.set_ylabel('Number of Oligosaccharides', 
                     fontsize=14)
        
        # 网格线设置
        ax.grid(alpha=0.3,
               color='gray',
               linestyle='--')
        
        # 自动调整y轴
        ax.set_ylim(0, max(values)*1.15 if values else 0)

    # 绘制子图
    create_plot(axes[0], single_dist, colors[0], 'Single Sugar Distribution')
    create_plot(axes[1], double_dist, colors[1], 'Double Sugar Distribution')
    create_plot(axes[2], triple_dist, colors[2], 'Triple Sugar Distribution')

    # 统一优化
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='x', which='both', length=3)  # 缩短刻度线

    plt.tight_layout()
    plt.show()

# 使用示例
if __name__ == "__main__":
    glycan_data = list(read_triplets('sugarlist.txt'))
    plot_sugar_distribution(glycan_data)