import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

datao = pd.read_csv('6.2.result-o.csv') # 普通算法
datas = pd.read_csv('6.2.result-s.csv') # 加速算法


length = 320
O1 = list(datao['M1'])[:length]
O2 = list(datao['M2'])[:length]
S1 = list(datas['M1'])[:length]
S2 = list(datas['M2'])[:length]

O1 = [np.log(x) for x in O1]
O2 = [np.log(x) for x in O2]
S1 = [np.log(x) for x in S1]
S2 = [np.log(x) for x in S2]

plt.rcParams['font.sans-serif']=['SimHei'] # 显示中文
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

plt.figure(figsize=(16,8))
plt.xlabel('iters')
plt.title(r"$\log(\rm{E}(\rm{error}_i))$")
plt.plot(O1,label='普通算法')
plt.plot(S1,label='加速算法')
plt.legend()
plt.show()

plt.figure(figsize=(16,8))
plt.xlabel('iters')
plt.title(r"$\log(\rm{E}(\rm{error}_i^2))$")
plt.plot(O2,label='普通算法')
plt.plot(S2,label='加速算法')
plt.legend()
plt.show()