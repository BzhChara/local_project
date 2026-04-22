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
t = np.linspace(0, 200, 100)  # 延长时间范围以观察长期动态

# 求解微分方程组
sol = odeint(model, Y0, t, args=(c, d, a, b, K, alpha), rtol=1e-6, atol=1e-6)

# 创建DataFrame
sol_df = pd.DataFrame(sol, columns=['H', 'P'])
sol_df['Time'] = t

# 绘制结果
plt.figure(figsize=(10, 6))
plt.plot(sol_df['Time'], sol_df['H'], label='Preys (H)')
plt.plot(sol_df['Time'], sol_df['P'], label='Predators (P)')
plt.xlabel('Time')
plt.ylabel('Numbers')
plt.title('Predator-Prey Model with Environmental Carrying Capacity and Predator Saturation')
plt.legend()
plt.grid(True)
plt.show()

# 保存结果到CSV文件
sol_df.to_csv('predator_prey_model_results.csv', index=False)