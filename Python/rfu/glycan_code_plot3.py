import re
import matplotlib.pyplot as plt

def count_branches(glycan_chain):
    # 使用正则表达式匹配所有括号内的内容，并统计括号对的数量
    branch_count = len(re.findall(r'\([^)]*\)', glycan_chain))
    return branch_count

def extract_branch_counts(filename):
    # 提取每条糖链的支链个数
    branch_counts = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        # Process each set of three lines as one glycan chain
        for i in range(0, len(lines), 3):
            # 检查是否有完整的三行
            if i + 2 >= len(lines):
                break
            single_line = lines[i].strip()
            if not single_line.startswith('single;'):
                continue  # 跳过无效的 `single` 行
            
            # 提取糖链结构
            parts = single_line.split(';')
            if len(parts) < 2:
                continue
            glycan_chain = parts[1]
            
            # 统计支链个数
            branch_count = count_branches(glycan_chain)
            branch_counts.append(branch_count)
    
    return branch_counts

def plot_branch_counts(branch_counts):
    # 设置图表样式与之前的示例一致
    plt.figure(figsize=(16, 8))  # 增大图表尺寸

    # 绘制折线图
    plt.plot(branch_counts, color='blue', label='Branch Count', marker='o', linestyle='-', linewidth=1, markersize=4, alpha=0.7)

    # 设置标题和坐标轴标签
    plt.title('Branch Counts of Glycan Chains', fontsize=16, fontweight='bold')
    plt.xlabel('Glycan Chain Index', fontsize=14)
    plt.ylabel('Number of Branches', fontsize=14)

    # 添加图例
    plt.legend(fontsize=12)

    # 调整横坐标显示
    plt.xticks(range(0, len(branch_counts), 10), rotation=45, fontsize=10)  # 每隔10个点显示一个横坐标

    # 调整网格线
    plt.grid(alpha=0.3, color='gray', linestyle='--')

    # 显示图表
    plt.tight_layout()
    plt.show()

# 示例用法
if __name__ == "__main__":
    input_file = 'sugarlist.txt'  # 替换为你的输入文件路径
    branch_counts = extract_branch_counts(input_file)
    plot_branch_counts(branch_counts)