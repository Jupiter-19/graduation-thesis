import os 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class Kubo():
    def __init__(self):
        self.X0, self.Y0 = 1,0
        self.b, self.sigma = 2,2
        self.h = 1e-2
        self.T = 80
        
        self.N = int(self.T / self.h)
        self.sqrt_h = np.sqrt(self.h)
        self.t = [self.h * i for i in range(self.N+1)]
        self.sqrt_I = np.sqrt(self.X0**2 + self.Y0**2)
        
    # 数据分布演化 num 次 
    def scheme(self,oldx,oldy,num):
        length = len(oldx)
        b,sigma,h = self.b , self.sigma, self.h
        
        newx,newy = [],[]
        for k in range(length):
            # 第 k 个点的模拟
            tx = [oldx[k] for _ in range(num+1)]
            ty = [oldy[k] for _ in range(num+1)]
            Brown = np.random.normal(0 , 1, self.N) * self.sqrt_h
            for i in range(num):
                # Taylor格式
                tx[i+1] = tx[i] - b*ty[i]*h - sigma*ty[i]*Brown[i] - sigma**2*tx[i]/2*Brown[i]**2
                tx[i+1] +=(Brown[i]**3-3*h*Brown[i])/6 * sigma**3*ty[i]
                tx[i+1] +=(sigma**4/4*tx[i]-sigma**2*b*ty[i]-b**2*tx[i])/2*h*h
                tx[i+1] -= Brown[i]*h*(sigma**3/2*ty[i]+b*sigma*tx[i])
                ty[i+1] = ty[i] + b*tx[i]*h + sigma*tx[i]*Brown[i] - sigma**2*ty[i]/2*Brown[i]**2
                ty[i+1] -=(Brown[i]**3-3*h*Brown[i])/6 * sigma**3*tx[i]
                ty[i+1] +=(sigma**4/4*ty[i]+sigma**2*b*tx[i]-b**2*ty[i])/2*h*h
                ty[i+1] += Brown[i]*h*(sigma**3/2*tx[i]-b*sigma*ty[i])
                # 投射操作
                tmp = np.sqrt(tx[i+1]**2 + ty[i+1]**2)
                tx[i+1] = tx[i+1]/tmp*self.sqrt_I
                ty[i+1] = ty[i+1]/tmp*self.sqrt_I
                
            newx.append(tx[-1])
            newy.append(ty[-1])

        return newx,newy
            
    def test(self):
        num = int(0.2 / self.h)
        
        for n in range(1,int(self.T/0.2)+1):
            time = str(n*2//10).zfill(2)+'.'+str(n*2%10)
            print(time)
            if 'T='+time+'.csv' in os.listdir('D:\data\Kubo\point'):
                data = pd.read_csv('D:\data\Kubo\point\T='+time+'.csv')
                x = list(data['x'])
                y = list(data['y'])
                continue
            
            x,y = self.scheme(x,y,num)
            out = pd.DataFrame({'x':x , 'y':y})
            out.to_csv('D:\data\Kubo\point\T='+time+'.csv')
            
            plt.hist(x,bins = 100)
            plt.savefig('D:\data\Kubo\histX\T='+time+'.png')
            plt.show()
            
            plt.hist(y,bins = 100)
            plt.savefig('D:\data\Kubo\histY\T='+time+'.png')
            plt.show()
            
            plt.figure(figsize=(8,8))
            plt.scatter(x[:500],y[:500])
            plt.show()
            #break
        
k = Kubo()
k.test()