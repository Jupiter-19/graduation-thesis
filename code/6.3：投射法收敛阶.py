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
        
        
        self.num = int(20000)
        self.N = int(self.T / self.h)
        
        self.sqrt_h = np.sqrt(self.h)
        self.t = [self.h * i for i in range(self.N+1)]


    def scheme(self, Ex,Ey,Mx,My,Tx,Ty,EPx,EPy,MPx,MPy ):
        b,sigma,h = self.b , self.sigma, self.h
        Brown = np.random.normal(0 , 1, self.num) * self.sqrt_h
        
        ex = [0 for i in range(self.num)]
        ey = [0 for i in range(self.num)]
        mx = [0 for i in range(self.num)]
        my = [0 for i in range(self.num)]
        tx = [0 for i in range(self.num)]
        ty = [0 for i in range(self.num)]
        
        epx= [0 for i in range(self.num)]
        epy= [0 for i in range(self.num)]
        mpx= [0 for i in range(self.num)]
        mpy= [0 for i in range(self.num)]
        
        
    
        
        for i in range(self.num):
            # Euler格式 Euler-Maruyama method
            ex[i] = Ex[i] - b*Ey[i]*h - sigma*Ey[i]*Brown[i] - sigma**2*Ex[i]/2*h
            ey[i] = Ey[i] + b*Ex[i]*h + sigma*Ex[i]*Brown[i] - sigma**2*Ey[i]/2*h
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
            # EulerP           
            epx[i] = EPx[i] - b*EPy[i]*h - sigma*EPy[i]*Brown[i] - sigma**2*EPx[i]/2*h
            epy[i] = EPy[i] + b*EPx[i]*h + sigma*EPx[i]*Brown[i] - sigma**2*EPy[i]/2*h            
            tmp = np.sqrt(epx[i]**2 + epy[i]**2)
            epx[i] = epx[i] / tmp
            epy[i] = epy[i] / tmp
            # MilsteinP
            mpx[i] = MPx[i] - b*MPy[i]*h - sigma*MPy[i]*Brown[i] - sigma**2*MPx[i]/2*Brown[i]**2
            mpy[i] = MPy[i] + b*MPx[i]*h + sigma*MPx[i]*Brown[i] - sigma**2*MPy[i]/2*Brown[i]**2
            tmp = np.sqrt(mpx[i]**2 + mpy[i]**2)
            mpx[i] = mpx[i] / tmp
            mpy[i] = mpy[i] / tmp
            
        return [ex,ey,mx,my,tx,ty,epx,epy,mpx,mpy]
            
    
    

    def test(self):
        Ex, Ey = [1]*self.num , [0]*self.num # Euler 格式
        Mx, My = [1]*self.num , [0]*self.num # Milstein 格式
        Tx, Ty = [1]*self.num , [0]*self.num # Taylor 格式
        EPx,EPy= [1]*self.num , [0]*self.num # Euler 格式
        MPx,MPy= [1]*self.num , [0]*self.num # Milstein 格式
        
        for idx in range(self.N):
            print(round(self.h*idx,4))
            output = self.scheme(Ex,Ey,Mx,My,Tx,Ty,EPx,EPy,MPx,MPy)
            Ex,Ey,Mx,My,Tx,Ty,EPx,EPy,MPx,MPy = output
        
        Ex, Ey = np.array(Ex) , np.array(Ey)
        Mx, My = np.array(Mx) , np.array(My)
        Tx, Ty = np.array(Tx) , np.array(Ty)
        
        EPx,EPy= np.array(EPx), np.array(EPy)
        MPx,MPy= np.array(MPx), np.array(MPy)
        
        out1 = [np.mean(np.abs(Ex-Tx)) , np.mean( (Ex-Tx)**2 )]
        out2 = [np.mean(np.abs(Mx-Tx)) , np.mean( (Mx-Tx)**2 )]
        out3 = [np.mean(np.abs(EPx-Tx)), np.mean( (EPx-Tx)**2 )]
        out4 = [np.mean(np.abs(MPx-Tx)), np.mean( (MPx-Tx)**2 )]
        
        output = out1+out2+out3+out4
        print(output)
        return output
 

# In[]
        
hlim = [1e-1 , 5e-2 , 2e-2 , 1e-2 , 5e-3 , 2e-3 , 1e-3 ]#, 5e-4]
headline = ['h','E1','E2','M1','M2','EP1','EP2','MP1','MP2']

out = open('6.3.result2.csv', 'a', newline='')
csv_write = csv.writer(out, dialect='excel')
csv_write.writerow(headline)
out.close()

for h_val in hlim:
    k = Kubo(h = h_val)
    output = k.test()
    
    out = open('6.3.result2.csv', 'a', newline='')
    csv_write = csv.writer(out, dialect='excel')
    csv_write.writerow([h_val] + output)
    
    out.close()


# In[]
data = pd.read_csv('C:/Users/Jupiter/Desktop/我的坚果云/研三上/.毕业设计/code/6.3.result2.csv')
h = [np.log(np.abs(item)) for item in data['h']]
E = [np.log(np.abs(item)) for item in data['E1']]
M = [np.log(np.abs(item)) for item in data['M1']]
EP= [np.log(np.abs(item)) for item in data['EP1']]
MP= [np.log(np.abs(item)) for item in data['MP1']]

plt.rcParams['font.sans-serif']=['SimHei'] # 显示中文
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

plt.figure(figsize=(16,8))
plt.xlabel('log(h)')
plt.title(r"$\log(\rm{E}(\rm{error}_i))$")
plt.plot(h,E ,label='Euler',     color='y',marker='*')
plt.plot(h,M ,label='Milstein',  color='r',marker='o')
plt.plot(h,EP,label='EulerP',    color='g',marker='v')
plt.plot(h,MP,label='MilsteinP', color='b',marker='s')
plt.legend() #显示图例
plt.show()


h = [np.log(np.abs(item)) for item in data['h']]
E = [np.log(np.abs(item)) for item in data['E2']]
M = [np.log(np.abs(item)) for item in data['M2']]
EP= [np.log(np.abs(item)) for item in data['EP2']]
MP= [np.log(np.abs(item)) for item in data['MP2']]

plt.rcParams['font.sans-serif']=['SimHei'] # 显示中文
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

plt.figure(figsize=(16,8))
plt.xlabel('log(h)')
plt.title(r"$\log(\rm{E}(\rm{error}_i^2))$")
plt.plot(h,E ,label='Euler',     color='y',marker='*')
plt.plot(h,M ,label='Milstein',  color='r',marker='o')
plt.plot(h,EP,label='EulerP',    color='g',marker='v')
plt.plot(h,MP,label='MilsteinP', color='b',marker='s')
plt.legend() #显示图例
plt.show()
