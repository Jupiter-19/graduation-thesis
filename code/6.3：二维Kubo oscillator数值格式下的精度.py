import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv

class Kubo():
    def __init__(self,h):
        self.X0, self.Y0 = 1 , 0
        self.b, self.sigma = 1 , 1
        self.h = h
        self.T = 0.5
        
        
        self.num = int(1e4)
        self.N = int(self.T / self.h)
        
        self.sqrt_h = np.sqrt(self.h)
        self.t = [self.h * i for i in range(self.N+1)]


    def scheme(self, E1x,E1y,E2x,E2y,Mx,My,Sx,Sy,Tx,Ty):
        b,sigma,h = self.b , self.sigma, self.h
        Brown = np.random.normal(0 , 1, self.num) * self.sqrt_h
        
        e1x = [0 for i in range(self.num)]
        e1y = [0 for i in range(self.num)]
        
        e2x = [0 for i in range(self.num)]
        e2y = [0 for i in range(self.num)]
        
        mx = [0 for i in range(self.num)]
        my = [0 for i in range(self.num)]
        
        sx = [0 for i in range(self.num)]
        sy = [0 for i in range(self.num)]
        
        tx = [0 for i in range(self.num)]
        ty = [0 for i in range(self.num)]
        
        for i in range(self.num):
            # Euler格式1 误差更大 
            e1x[i] = E1x[i] - b*E1y[i]*h - sigma*E1y[i]*Brown[i]
            e1y[i] = E1y[i] + b*E1x[i]*h + sigma*E1x[i]*Brown[i]
            # Euler格式2 更正确 (Euler-Maruyama method)
            e2x[i] = E2x[i] - b*E2y[i]*h - sigma*E2y[i]*Brown[i] - sigma**2*E1x[i]/2*h
            e2y[i] = E2y[i] + b*E2x[i]*h + sigma*E2x[i]*Brown[i] - sigma**2*E2y[i]/2*h
            # Milstein格式
            mx[i] = Mx[i] - b*My[i]*h - sigma*My[i]*Brown[i] - sigma**2*Mx[i]/2*Brown[i]**2
            my[i] = My[i] + b*Mx[i]*h + sigma*Mx[i]*Brown[i] - sigma**2*My[i]/2*Brown[i]**2
            # Taylor格式
            tx[i] = Tx[i] - b*Ty[i]*h - sigma*Ty[i]*Brown[i] - sigma**2*Tx[i]/2*Brown[i]**2
            tx[i] +=(Brown[i]**3-3*h*Brown[i])/6 * sigma**3*Ty[i]
            tx[i] +=(sigma**4/4*Tx[i]-sigma**2*b*Ty[i]-b**2*Tx[i])/2*h*h
            tx[i] -= Brown[i]*h*(sigma**3/2*Ty[i]+b*sigma*Tx[i])
            ty[i] = Ty[i] + b*Tx[i]*h + sigma*Tx[i]*Brown[i] - sigma**2*Ty[i]/2*Brown[i]**2
            ty[i] -=(Brown[i]**3-3*h*Brown[i])/6 * sigma**3*Tx[i]
            ty[i] +=(sigma**4/4*Ty[i]+sigma**2*b*Tx[i]-b**2*Ty[i])/2*h*h
            ty[i] += Brown[i]*h*(sigma**3/2*Tx[i]-b*sigma*Ty[i])
            # 保守恒量格式
            tmp = b*h/2 + sigma*Brown[i]/2
            a = Sx[i] - tmp*Sy[i]
            b = tmp*Sx[i] + Sy[i]
            sx[i] = (a-tmp*b)/(1+tmp**2)
            sy[i] = (tmp*a+b)/(1+tmp**2)

        return [e1x,e2x,mx,sx,tx, e1y,e2y,my,sy,ty]
            
    
    

    def test(self):
        E1x,E1y = [1]*self.num , [0]*self.num # 直接差分
        E2x,E2y = [1]*self.num , [0]*self.num # Euler 格式
        Mx, My  = [1]*self.num , [0]*self.num # Milstein 格式
        Sx, Sy  = [1]*self.num , [0]*self.num # 保守恒量格式
        Tx, Ty  = [1]*self.num , [0]*self.num # Taylor 格式
        
        for idx in range(self.N):
            print(round(self.h*idx,4))
            output = self.scheme(E1x,E1y,E2x,E2y,Mx,My,Sx,Sy,Tx,Ty)
            E1x,E2x,Mx,Sx,Tx, E1y,E2y,My,Sy,Ty = output
            
        
        I1 = np.array([ E1x[i]**2+E1y[i]**2 for i in range(self.num)])
        I2 = np.array([ E2x[i]**2+E2y[i]**2 for i in range(self.num)])
        I3 = np.array([ Mx[i]**2+My[i]**2 for i in range(self.num)])
        I4 = np.array([ Sx[i]**2+Sy[i]**2 for i in range(self.num)])
        I5 = np.array([ Tx[i]**2+Ty[i]**2 for i in range(self.num)])
        
        error1 = np.abs(np.array([ E1x[i] - Tx[i] for i in range(self.num) ]))
        error2 = np.abs(np.array([ E2x[i] - Tx[i] for i in range(self.num) ]))
        error3 = np.abs(np.array([ Mx[i]  - Tx[i] for i in range(self.num) ]))
        error4 = np.abs(np.array([ Sx[i]  - Tx[i] for i in range(self.num) ]))
        
        out1 = [np.mean(error1), np.mean(error1**2), np.mean(I1), np.std(I1)]
        out2 = [np.mean(error2), np.mean(error2**2), np.mean(I2), np.std(I2)]
        out3 = [np.mean(error3), np.mean(error3**2), np.mean(I3), np.std(I3)]
        out4 = [np.mean(error4), np.mean(error4**2), np.mean(I4), np.std(I4)]
        out5 = [np.mean(I5), np.std(I5)]
        
        
        return  out1 + out2 + out3 + out4 + out5
 
hlim = [1e-1 , 5e-2 , 2e-2 , 1e-2 , 5e-3 , 2e-3 , 1e-3 , 5e-4]           


# In[]

out = pd.DataFrame(columns=('h','E1L1','E1L2','E1IL1','E1IL2',
                                'E2L1','E2L2','E2IL1','E2IL2',
                                'ML1', 'ML2', 'MIL1', 'MIL2',
                                'SL1', 'SL2', 'SIL1', 'SIL2',
                                'TIL1', 'TIL2'))


headline = ['h','E1L1','E1L2','E1IL1','E1IL2','E2L1','E2L2','E2IL1','E2IL2',
            'ML1', 'ML2', 'MIL1', 'MIL2','SL1', 'SL2', 'SIL1', 'SIL2','TIL1', 'TIL2']

out = open('6.3.result1.csv', 'a', newline='')
csv_write = csv.writer(out, dialect='excel')
csv_write.writerow(headline)
out.close()

for h_val in hlim:
    k = Kubo(h = h_val)
    output = k.test()
    
    out = open('6.3.result1.csv', 'a', newline='')
    csv_write = csv.writer(out, dialect='excel')
    csv_write.writerow([h_val] + output)
    
    out.close()

# In[] 结果呈现
data = pd.read_csv('C:/Users/Jupiter/Desktop/我的坚果云/研三上/.毕业设计/code/6.3.result1.csv')
h = [np.log(np.abs(item)) for item in data['h']]
E1= [np.log(np.abs(item)) for item in data['E1L1']]
E2= [np.log(np.abs(item)) for item in data['E2L1']]
M = [np.log(np.abs(item)) for item in data['ML1']]
S = [np.log(np.abs(item)) for item in data['SL1']]

plt.rcParams['font.sans-serif']=['SimHei'] # 显示中文
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

plt.figure(figsize=(16,8))
plt.xlabel('log(h)')
plt.title(r"$\log(\rm{E}(\rm{error}_i))$")
plt.plot(h,E1,label='直接差分',     color='y',marker='*')
plt.plot(h,E2,label='Euler',       color='r',marker='o')
plt.plot(h,M, label='Milstein',    color='g',marker='v')
plt.plot(h,S, label='保守恒量格式',  color='b',marker='s')
plt.legend() #显示图例
plt.show()


h = [np.log(np.abs(item)) for item in data['h']]
E1= [np.log(np.abs(item)) for item in data['E1L2']]
E2= [np.log(np.abs(item)) for item in data['E2L2']]
M = [np.log(np.abs(item)) for item in data['ML2']]
S = [np.log(np.abs(item)) for item in data['SL2']]

plt.rcParams['font.sans-serif']=['SimHei'] # 显示中文
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

plt.figure(figsize=(16,8))
plt.xlabel('log(h)')
plt.title(r"$\log(\rm{E}(\rm{error}_i^2))$")
plt.plot(h,E1,label='直接差分',     color='y',marker='*')
plt.plot(h,E2,label='Euler',       color='r',marker='o')
plt.plot(h,M, label='Milstein',    color='g',marker='v')
plt.plot(h,S, label='保守恒量格式',  color='b',marker='s')
plt.legend() #显示图例
plt.show()


h = [np.log(np.abs(item)) for item in data['h']]
E1= [np.abs(item) for item in data['E1IL1']]
E2= [np.abs(item) for item in data['E2IL1']]
M = [np.abs(item) for item in data['MIL1']]
S = [np.abs(item) for item in data['SIL1']]
T = [np.abs(item) for item in data['TIL1']]

plt.rcParams['font.sans-serif']=['SimHei'] # 显示中文
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

plt.figure(figsize=(16,8))
plt.xlabel('log(h)')
plt.title('守恒量的期望')
plt.plot(h,E1,label='直接差分',     color='y',marker='*')
plt.plot(h,E2,label='Euler',       color='r',marker='o')
plt.plot(h,M, label='Milstein',    color='g',marker='v')
plt.plot(h,S, label='保守恒量格式', color='b',marker='s')
plt.plot(h,T, label='Taylor',      color='gray',marker='^')
plt.legend() #显示图例
plt.show()


h = [np.log(np.abs(item)) for item in data['h']]
E1= [np.abs(item) for item in data['E1IL2']]
E2= [np.abs(item) for item in data['E2IL2']]
M = [np.abs(item) for item in data['MIL2']]
S = [np.abs(item) for item in data['SIL2']]
T = [np.abs(item) for item in data['TIL2']]

plt.rcParams['font.sans-serif']=['SimHei'] # 显示中文
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

plt.figure(figsize=(16,8))
plt.xlabel('log(h)')
plt.title('守恒量的方差')
plt.plot(h,E1,label='直接差分',     color='y',marker='*')
plt.plot(h,E2,label='Euler',       color='r',marker='o')
plt.plot(h,M, label='Milstein',    color='g',marker='v')
plt.plot(h,S, label='保守恒量格式', color='b',marker='s')
plt.plot(h,T, label='Taylor',      color='gray',marker='^')
plt.legend() #显示图例
plt.show()

