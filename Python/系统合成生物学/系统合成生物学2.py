import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# 定义微分方程组，考虑猎物的环境承载力和捕食者的饱和效应
def model(Y, t, c, d, a, b, K, alpha):
    H, P = Y
    # 猎物的增长考虑环境承载力和捕食者的影响
    dH_dt = c * H * (1 - H / K) - d * H * P
    # 捕食者的增长包含饱和效应
    dP_dt = a * (H / (1 + alpha * H)) * P - b * P
    return [dH_dt, dP_dt]

# 参数设定
K = 100  # 环境承载力
c = 0.5  # 猎物的出生率
d = 0.1  # 捕食者与猎物之间的相互作用强度（捕食率）
a = 0.05 # 捕食者与猎物之间的相互作用强度（捕食效率）
b = 0.3  # 捕食者的死亡率
alpha = 0.01  # 捕食者饱和常数

# 初始条件
H0 = 50  # 初始猎物数量
P0 = 10  # 初始捕食者数量
Y0 = [H0, P0]

# 时间范围
t = np.linspace(0, 400, 100)  # 延长时间范围以观察长期动态

# 猎物出生率的变化范围和步长
c_range = np.arange(0.1, 1.1, 0.1)  # 从0.1到1.0，步长为0.1

# 创建一个空的列表来存储所有结果
all_results = []

for c in c_range:
    # 求解微分方程组
    sol = odeint(model, Y0, t, args=(c, d, a, b, K, alpha), rtol=1e-6, atol=1e-6)
    # 将解转换为DataFrame，并添加出生率列
    sol_df = pd.DataFrame(sol, columns=['H', 'P'])
    sol_df['Time'] = t
    sol_df['Birth Rate (c)'] = c
    # 将结果存储到列表中
    all_results.append(sol_df)

# 使用concat合并所有结果为一个DataFrame
df = pd.concat(all_results, ignore_index=True)

# 绘制结果
plt.figure(figsize=(10, 6))
for c_value in c_range:
    subset = df[(df['Birth Rate (c)'] == c_value)]
    plt.plot(subset['Time'], subset['H'], label=f'Birth Rate: {c_value}')

plt.xlabel('Time')
plt.ylabel('Numbers')
plt.title('Predator-Prey Model with Varying Birth Rate (c)')
plt.legend()
plt.grid(True)
plt.show()

# 保存结果到CSV文件
df.to_csv('predator_prey_model_results.csv', index=False)