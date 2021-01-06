import numpy as np
import matplotlib.pyplot as plt

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
        
        
        I11,I21,I31,I41,I51 = [1],[1],[1],[1],[1]
        I12,I22,I32,I42,I52 = [0],[0],[0],[0],[0]
        
        for idx in range(self.N):
            print(round(self.h*idx,4))
            output = self.scheme(E1x,E1y,E2x,E2y,Mx,My,Sx,Sy,Tx,Ty)
            E1x,E2x,Mx,Sx,Tx, E1y,E2y,My,Sy,Ty = output
            
        
            I1 = np.array([ E1x[i]**2+E1y[i]**2 for i in range(self.num)])
            I2 = np.array([ E2x[i]**2+E2y[i]**2 for i in range(self.num)])
            I3 = np.array([ Mx[i]**2+My[i]**2 for i in range(self.num)])
            I4 = np.array([ Sx[i]**2+Sy[i]**2 for i in range(self.num)])
            I5 = np.array([ Tx[i]**2+Ty[i]**2 for i in range(self.num)])
            
            I11.append(np.mean(I1))
            I21.append(np.mean(I2))
            I31.append(np.mean(I3))
            I41.append(np.mean(I4))
            I51.append(np.mean(I5))
            
            I12.append(np.std(I1))
            I22.append(np.std(I2))
            I32.append(np.std(I3))
            I42.append(np.std(I4))
            I52.append(np.std(I5))
        
        return self.t,I11,I21,I31,I41,I51,I12,I22,I32,I42,I52 
        
 

k = Kubo(h = 5e-2)
t,I11,I21,I31,I41,I51,I12,I22,I32,I42,I52 = k.test()
# In[]

plt.rcParams['font.sans-serif']=['SimHei'] # 显示中文
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

plt.figure(figsize=(16,8))
plt.xlabel('t')
plt.title('守恒量的均值随时间的变化(h = 0.05)')
plt.plot(t,I11,label='直接差分',     color='y')
plt.plot(t,I21,label='Euler',       color='r')
plt.plot(t,I31,label='Milstein',    color='g')
plt.plot(t,I41,label='保守恒量格式', color='b')
plt.plot(t,I51,label='Taylor',      color='gray')
plt.legend() #显示图例
plt.show()



plt.rcParams['font.sans-serif']=['SimHei'] # 显示中文 
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

plt.figure(figsize=(16,8))
plt.xlabel('t')
plt.title('守恒量的方差随时间的变化(h = 0.05)')
plt.plot(t,I12,label='直接差分',    color='y')
plt.plot(t,I22,label='Euler',       color='r')
plt.plot(t,I32,label='Milstein',    color='g')
plt.plot(t,I42,label='保守恒量格式', color='b')
plt.plot(t,I52,label='Taylor',      color='gray')
plt.legend() #显示图例
plt.show()