import re
import matplotlib.pyplot as plt

def read_triplets(file_path):
    """生成器函数，每次读取三行作为一组"""
    with open(file_path, 'r') as file:
        lines = []
        for line in file:
            lines.append(line.strip())
            if len(lines) == 3:
                yield lines
                lines = []

def extract_total_sugars(line):
    """提取糖类总数"""
    counts = line.split(';')[-1]
    numbers = re.findall(r'(\d+)', counts)
    return sum(map(int, numbers))

def plot_sugar_counts(glycan_data):
    """绘制糖类数量折线图"""
    single_counts = []
    double_counts = []
    trip_counts = []
    
    for triplet in glycan_data:
        single_line = triplet[0]
        double_line = triplet[1]
        trip_line = triplet[2]
        
        single_count = extract_total_sugars(single_line)
        double_count = extract_total_sugars(double_line)
        trip_count = extract_total_sugars(trip_line)
        
        single_counts.append(single_count)
        double_counts.append(double_count)
        trip_counts.append(trip_count)
    
    # 绘图
    plt.figure(figsize=(16, 8))  # 增大图表尺寸
    
    # 绘制单糖、二糖、三糖的折线图
    plt.plot(single_counts, color='blue', label='Single Sugars', marker='o', linestyle='-', linewidth=1, markersize=4, alpha=0.7)
    plt.plot(double_counts, color='red', label='Double Sugars', marker='s', linestyle='-', linewidth=1, markersize=4, alpha=0.7)
    plt.plot(trip_counts, color='green', label='Triple Sugars', marker='^', linestyle='-', linewidth=1, markersize=4, alpha=0.7)
    
    # 设置标题和坐标轴标签
    plt.title('Sugar Counts in Glycans', fontsize=16, fontweight='bold')
    plt.xlabel('Glycan Index', fontsize=14)
    plt.ylabel('Number of Sugars', fontsize=14)
    
    # 添加图例
    plt.legend(fontsize=12)
    
    # 调整横坐标显示
    plt.xticks(range(0, len(single_counts), 10), rotation=45, fontsize=10)  # 每隔10个点显示一个横坐标
    
    # 调整网格线
    plt.grid(alpha=0.3, color='gray', linestyle='--')  # 调整网格线透明度和颜色
    
    # 显示图表
    plt.tight_layout()
    plt.show()

# 示例使用
if __name__ == "__main__":
    # 假设数据存储在名为 glycan_data.txt 的文件中
    glycan_data = list(read_triplets('sugarlist.txt'))
    
    # 绘制折线图
    plot_sugar_counts(glycan_data)