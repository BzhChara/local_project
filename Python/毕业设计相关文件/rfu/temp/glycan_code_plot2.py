import re
import matplotlib.pyplot as plt

def extract_main_chain_lengths(filename):
    # Extract main chain lengths from the glycan data file
    chain_lengths = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        # Process each set of three lines as one glycan chain
        for i in range(0, len(lines), 3):
            # Check if we have a complete group of three lines
            if i + 2 >= len(lines):
                break
            single_line = lines[i].strip()
            if not single_line.startswith('single;'):
                continue  # Skip if not a valid single line
            
            # Extract the glycan structure from the single line
            parts = single_line.split(';')
            if len(parts) < 2:
                continue
            glycan_chain = parts[1]
            
            # Remove bracketed sections (branch chains) to get the main chain
            main_chain = re.sub(r'\([^)]*\)', '', glycan_chain)
            
            # Extract sugar names using regex to find all sequences of letters
            sugar_names = re.findall(r'[A-Za-z]+', main_chain)
            
            # Count the number of sugar units on the main chain
            chain_length = len(sugar_names)
            chain_lengths.append(chain_length)
    
    return chain_lengths

def plot_main_chain_lengths(chain_lengths):
    plt.figure(figsize=(16, 8))  # 增大图表尺寸

    # 绘制主链长度折线图
    plt.plot(chain_lengths, color='blue', label='Main Chain Length', marker='o', linestyle='-', linewidth=1, markersize=4, alpha=0.7)

    # 设置标题和坐标轴标签
    plt.title('Main Chain Lengths of Glycan Chains', fontsize=16, fontweight='bold')
    plt.xlabel('Glycan Chain Index', fontsize=14)
    plt.ylabel('Main Chain Length', fontsize=14)

    # 添加图例
    plt.legend(fontsize=12)

    # 调整横坐标显示
    plt.xticks(range(0, len(chain_lengths), 10), rotation=45, fontsize=10)  # 每隔10个点显示一个横坐标

    # 调整网格线
    plt.grid(alpha=0.3, color='gray', linestyle='--')  # 调整网格线透明度和颜色

    # 显示图表
    plt.tight_layout()
    plt.show()

# 示例用法
if __name__ == "__main__":
    input_file = 'sugarlist.txt'  # 替换为你的输入文件路径
    chain_lengths = extract_main_chain_lengths(input_file)
    plot_main_chain_lengths(chain_lengths)