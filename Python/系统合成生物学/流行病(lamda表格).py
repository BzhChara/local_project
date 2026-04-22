# SEIR 模型，常微分方程组
from scipy.integrate import odeint  # 导入 scipy.integrate 模块
import numpy as np  # 导入 numpy包
import matplotlib.pyplot as plt  # 导入 matplotlib包
import pandas as pd 

def dySEIR(y, t, lamda, delta, mu):  # SEIR 模型，导数函数
    s, e, i = y  
    ds_dt = -lamda*s*i  
    de_dt = lamda*s*i - delta*e  
    di_dt = delta*e - mu*i  
    return np.array([ds_dt,de_dt,di_dt])

# 设置模型参数
number = 1e5  # 总人数
delta = 0.1  # 日发病率，每天发病成为患病者的潜伏者占潜伏者总数的比例
mu = 0.06  # 日治愈率, 每天治愈的患病者人数占患病者总数的比例

tEnd = 300  # 预测日期长度
t = np.arange(0.0,tEnd,1)    
i0 = 1e-3  # 患病者比例的初值
e0 = 2e-3  # 潜伏者比例的初值
s0 = 0.997  # 易感者比例的初值

lamda_values = [0.12, 0.25, 0.5, 1.0, 2.0]  # 日接触率, 患病者每天有效接触的易感者的平均人数

colors = ['m', 'c', 'y', 'k', 'r']  # 对应的颜色列表

# 初始化一个空列表来存储结果
results_list = []

# 输出绘图
for lamda in lamda_values:

    sigma = lamda / mu  # 传染期接触数
    fsig = 1-1/sigma

    Y0 = (s0, e0, i0)  # 微分方程组的初值
    ySEIR = odeint(dySEIR, Y0, t, args=(lamda, delta, mu))
    
    # 从结果中计算所需的指标
    i_max = np.max(ySEIR[:, 2])  # 患病者比例峰值
    s_final = ySEIR[-1, 0]  # 易感者比例最终稳定值
    t_peak_index = np.argmax(ySEIR[:, 2])  # 疫情达到峰值的索引
    t_peak = t[t_peak_index]  # 疫情达到峰值的时间

    # 将结果以字典形式添加到列表中
    results_list.append({
        'λ': lamda,
        'i_max': i_max,
        's_final': s_final,
        't_peak': t_peak
    })

    # 使用列表创建DataFrame
    results_df = pd.DataFrame(results_list)

    col=colors.pop(0)

    # 绘制曲线
    plt.plot(t, ySEIR[:,0], '--', color=col, label=f'λ={lamda:0%}'+'s(t)')
    plt.plot(t, ySEIR[:,2], '-', color=col, label=f'λ={lamda:.0%}'+'i(t)')

# 输出表格到控制台
print(results_df)

# 添加图表元素
plt.title('lmpact of λ on i(t),s(t) in SElR model')
plt.xlabel('Time (days)')
plt.ylabel('Proportion of Population')
plt.legend(loc='best')  # 显示图例
plt.grid(True)  # 显示网格
plt.tight_layout()  # 自动调整子图布局以防止重叠
plt.show()