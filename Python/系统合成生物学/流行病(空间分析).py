# SEIR 模型，常微分方程组 相空间分析: e(t)~i(t)
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
delta = 0.1  # 日发病率，每天发病成为患病者的潜伏者占潜伏者总数的比例
mu = 0.1  # 日治愈率, 每天治愈的患病者人数占患病者总数的比例
sigma = lamda / mu  # 传染期接触数
tEnd = 500  # 预测日期长度
t = np.arange(0.0, tEnd, 1)  

e0List = np.arange(0.01,0.4,0.05)  # (start,stop,step)
for e0 in e0List:
    # odeint 数值解，求解微分方程初值问题
    i0 = 0  # 潜伏者比例的初值
    s0 = 1 - i0 - e0  # 易感者比例的初值
    ySEIR = odeint(dySEIR, (s0,e0,i0), t, args=(lamda,delta,mu))  # SEIR 模型
    plt.plot(ySEIR[:,0], ySEIR[:,1])  # plt.plot(ySEIR[:,1], ySEIR[:,2])
    print("lamda={}\tdelta={}\mu={}\tsigma={}\ti0={}\te0={}".format(lamda,delta,mu,lamda/mu,i0,e0))

# 输出绘图
plt.title("Phase trajectory of SEIR models: s(t)~e(t)")
# plt.title("Phase trajectory of SEIR models: e(t)~i(t)")
plt.axis([0, 1.0, 0, 1.0])
# plt.axis([0, 0.4, 0, 0.4])
plt.plot([0,1.0],[1.0,0],'b-')  #[x1,x2][y1,y2]
plt.plot([0,1.0],[0.8,0],'y--')  #[x1,x2][y1,y2]
# plt.plot([0,0.4],[0,0.35],'y--')  
# plt.plot([0,0.4],[0,0.18],'y--')
plt.text(0.6,0.8,r"$\lambda=0.3, \delta=0.1, \mu=0.1$",color='black')
# plt.text(0.02,0.36,r"$\lambda=0.3, \delta=0.1, \mu=0.1$",color='black')
plt.xlabel('s(t)')
# plt.xlabel('e(t)')
plt.ylabel('e(t)')
# plt.ylabel('i(t)')
plt.show()