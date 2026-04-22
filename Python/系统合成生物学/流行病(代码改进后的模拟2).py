# SEIR2 模型，考虑潜伏期具有传染性
from scipy.integrate import odeint  # 导入 scipy.integrate 模块
import numpy as np  # 导入 numpy包
import matplotlib.pyplot as plt  # 导入 matplotlib包

def dySEIR(y, t, lamda, delta, mu):  # SEIR 模型，导数函数
    s, e, i = y
    ds_dt = - lamda*s*i  
    de_dt = lamda*s*i - delta*e  
    di_dt = delta*e - mu*i  
    return np.array([ds_dt,de_dt,di_dt])

def dySEIR2(y, t, lamda, lam2, delta, mu):  # SEIR2 模型，导数函数
    s, e, i = y
    ds_dt = - lamda*s*i - lam2*s*e 
    de_dt = lamda*s*i + lam2*s*e - delta*e  
    di_dt = delta*e - mu*i 
    return np.array([ds_dt,de_dt,di_dt])

# 设置模型参数
number = 1e5  # 总人数

delta = 0.05  # 日发病率，每天发病成为患病者的潜伏者占潜伏者总数的比例
mu = 0.05  # 日治愈率, 每天治愈的患病者人数占患病者总数的比例

tEnd = 200  # 预测日期长度
t = np.arange(0.0, tEnd, 1)  # (start,stop,step)
i0 = 1e-3  # 患病者比例的初值
e0 = 0  # 潜伏者比例的初值
s0 = 1-i0  # 易感者比例的初值
Y0 = (s0, e0, i0)  # 微分方程组的初值

lamda = [1.0, 0.2, 0.2]  # 日接触率, 患病者每天有效接触的易感者的平均人数
lam2 = [0.25, 0.25, 0.1]  # 日接触率2, 潜伏者每天有效接触的易感者的平均人数

linestyle_values = ['-', '--', ':']

for num in range(3):
    sigma = lamda[num] / mu  # 传染期接触数
    fsig = 1-1/sigma

    line = linestyle_values.pop(0)

    # odeint 数值解，求解微分方程初值问题
    ySEIR2 = odeint(dySEIR2, Y0, t, args=(lamda[num],lam2[num],delta,mu))  # SEIR2 模型

    # 输出绘图
    plt.plot(t, ySEIR2[:,0], 'g', label='s(t)-iSEIR', linestyle = line)  # 易感者比例
    plt.plot(t, ySEIR2[:,1], 'b', label='e(t)-iSEIR', linestyle = line)  # 潜伏者比例
    plt.plot(t, ySEIR2[:,2], 'm', label='i(t)-iSEIR', linestyle = line)  # 患病者比例

# 添加图表元素
plt.title("Comparison among different prevention methods")
plt.xlabel('t')
plt.axis([0, tEnd, -0.1, 1.1])
plt.legend(loc='best')  # 显示图例
plt.show()