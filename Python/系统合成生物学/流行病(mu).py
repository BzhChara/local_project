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
lamda = 0.25  # 日接触率, 患病者每天有效接触的易感者的平均人数
delta = 0.1 # 日发病率，每天发病成为患病者的潜伏者占潜伏者总数的比例
  
tEnd = 300  # 预测日期长度
t = np.arange(0.0,tEnd,1)  
i0 = 1e-3  # 患病者比例的初值
e0 = 2e-3  # 潜伏者比例的初值
s0 = 0.997  # 易感者比例的初值
Y0 = (s0, e0, i0)  # 微分方程组的初值

mu_values = [0.02, 0.05, 0.1, 0.2, 0.4] # 日治愈率, 每天治愈的患病者人数占患病者总数的比例

colors = ['m', 'c', 'y', 'k', 'r']  # 对应的颜色列表

# 输出绘图
for mu in mu_values:
    sigma = lamda / mu  # 传染期接触数
    fsig = 1-1/sigma

    Y0 = (s0, e0, i0)  # 微分方程组的初值
    ySEIR = odeint(dySEIR, Y0, t, args=(lamda, delta, mu))
    
    col=colors.pop(0)

    # 绘制曲线
    plt.plot(t, ySEIR[:,0], '--', color=col, label=f'mu={mu:.0%}'+'s(t)')
    plt.plot(t, ySEIR[:,2], '-', color=col, label=f'mu={mu:.0%}'+'i(t)')

# 添加图表元素
plt.title('lmpact of μ on i(t),s(t) in SElR model')
plt.xlabel('t')
plt.ylabel('Proportion of Population')
plt.legend(loc='best')  # 显示图例
plt.grid(True)  # 显示网格
plt.tight_layout()  # 自动调整子图布局以防止重叠
plt.show()