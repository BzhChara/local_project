# SEIR 模型，常微分方程组
from scipy.integrate import odeint  # 导入 scipy.integrate 模块
import numpy as np  # 导入 numpy包
import matplotlib.pyplot as plt  # 导入 matplotlib包

def dySEIR(y, t, lamda, delta, mu):  # SEIR 模型，导数函数
    s, e, i = y  
    ds_dt = -lamda*s*i  
    de_dt = lamda*s*i - delta*e  
    di_dt = delta*e - mu*i  
    return np.array([ds_dt,de_dt,di_dt])

# 设置模型参数
number = 1e5  # 总人数
lamda = 0.3  # 日接触率, 患病者每天有效接触的易感者的平均人数
delta = 0.03  # 日发病率，每天发病成为患病者的潜伏者占潜伏者总数的比例
mu = 0.06  # 日治愈率, 每天治愈的患病者人数占患病者总数的比例
sigma = lamda / mu  # 传染期接触数
fsig = 1-1/sigma
tEnd = 300  # 预测日期长度
t = np.arange(0.0,tEnd,1)    

i0_values = [0.1, 0.03, 0.01, 0.001, 0.0001] # 患病者比例的初值

colors = ['m', 'c', 'y', 'k', 'r']  # 对应的颜色列表

# 输出绘图
for i0 in i0_values:
    s0 = 1-i0
    e0 = 2*i0  
    Y0 = (s0, e0, i0)  # 微分方程组的初值
    ySEIR = odeint(dySEIR, Y0, t, args=(lamda, delta, mu))
    
    col=colors.pop(0)

    # 绘制曲线
    plt.plot(t, ySEIR[:,0], '--', color=col, label=f'i0={i0:.0%}'+'s(t)')
    plt.plot(t, ySEIR[:,2], '-', color=col, label=f'i0={i0:.0%}'+'i(t)')

# 添加图表元素
plt.title('lmpact of i0 on i(t),s(t) in SElR model')
plt.xlabel('t')
plt.ylabel('Proportion of Population')
plt.legend(loc='best')  # 显示图例
plt.grid(True)  # 显示网格
plt.tight_layout()  # 自动调整子图布局以防止重叠
plt.show()