import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root,fsolve
from mpl_toolkits.mplot3d import Axes3D 

class Lotka():
    def __init__(self,h,N):
        self.h = h
        self.N = N
        self.I1 = 4
        self.I2 = 2
        
        self.num = 5000 # 样本轨道数

    def loss(self,x):
        return np.array([x[0]-x[3]-x[4]*x[1]*x[2]-self.a,
                         x[1]-x[3]-x[4]*x[0]*x[2]-self.b,
                         x[2]-x[3]-x[4]*x[0]*x[1]-self.c,
                         x[0]+x[1]+x[2]-self.I1,
                         x[0]*x[1]*x[2]-self.I2])
    
    
    def scheme(self, Ex,Ey,Ez , Mx,My,Mz , EPx,EPy,EPz , MPx,MPy,MPz):
        h = self.h
        Brown = np.random.normal(0 , 1, self.num) * np.sqrt(h)
        
        ex = [0 for i in range(self.num)]
        ey = [0 for i in range(self.num)]
        ez = [0 for i in range(self.num)]
        mx = [0 for i in range(self.num)]
        my = [0 for i in range(self.num)]
        mz = [0 for i in range(self.num)]
        
        epx= [0 for i in range(self.num)]
        epy= [0 for i in range(self.num)]
        epz= [0 for i in range(self.num)]
        mpx= [0 for i in range(self.num)]
        mpy= [0 for i in range(self.num)]
        mpz= [0 for i in range(self.num)]
        
        
    
        for i in range(self.num):
            # Euler格式 Euler-Maruyama method
            # Milstein格式
            mx[i] = Mx[i] + Mx[i]*(Mz[i]-My[i])*(h+Brown[i])
            mx[i]+= Brown[i]**2/2 * Mx[i] * (My[i]**2+Mz[i]**2-Mx[i]*(My[i]+Mz[i]))
            
            my[i] = My[i] + My[i]*(Mx[i]-Mz[i])*(h+Brown[i])
            my[i]+= Brown[i]**2/2 * My[i] * (Mx[i]**2+Mz[i]**2-My[i]*(Mx[i]+Mz[i]))
            
            mz[i] = Mz[i] + Mz[i]*(My[i]-Mx[i])*(h+Brown[i])
            mz[i]+= Brown[i]**2/2 * Mz[i] * (Mx[i]**2+My[i]**2-Mz[i]*(Mx[i]+My[i]))            
            # EulerP           
            # MilsteinP
            mpx[i] = MPx[i] + MPx[i]*(MPz[i]-MPy[i])*(h+Brown[i])
            mpx[i]+= Brown[i]**2/2 * MPx[i] * (MPy[i]**2+MPz[i]**2-MPx[i]*(MPy[i]+MPz[i]))
            
            mpy[i] = MPy[i] + MPy[i]*(MPx[i]-MPz[i])*(h+Brown[i])
            mpy[i]+= Brown[i]**2/2 * MPy[i] * (MPx[i]**2+MPz[i]**2-MPy[i]*(MPx[i]+MPz[i]))
            
            mpz[i] = MPz[i] + MPz[i]*(MPy[i]-MPx[i])*(h+Brown[i])
            mpz[i]+= Brown[i]**2/2 * MPz[i] * (MPx[i]**2+MPy[i]**2-MPz[i]*(MPx[i]+MPy[i])) 
            
            self.a, self.b, self.c = mpx[i],mpy[i],mpz[i]
            sol_root = root(self.loss,[self.a, self.b, self.c ,0 ,0]).x
            mpx[i],mpy[i],mpz[i] = sol_root[0] , sol_root[1] , sol_root[2]

        return [ex,ey,ez , mx,my,mz , epx,epy,epz , mpx,mpy,mpz]
            
    
    

    def test(self):
        Ex, Ey, Ez = [2]*self.num , [1]*self.num , [1]*self.num # Euler 格式
        Mx, My, Mz = [2]*self.num , [1]*self.num , [1]*self.num # Milstein 格式
        EPx,EPy,EPz= [2]*self.num , [1]*self.num , [1]*self.num # EulerP 格式
        MPx,MPy,MPz= [2]*self.num , [1]*self.num , [1]*self.num # MilsteinP 格式
        
        for idx in range(self.N):
            print(round(self.h*idx,4))

            output = self.scheme(Ex,Ey,Ez , Mx,My,Mz , EPx,EPy,EPz , MPx,MPy,MPz)
            Ex,Ey,Ez,Mx,My,Mz,EPx,EPy,EPz,MPx,MPy,MPz = output

        fig = plt.figure( figsize=(12,8) )
        ax = Axes3D(fig)
        ax.scatter(MPx,MPy,MPz, c='r', label='MilsteinP', s=1)     
        ax.scatter(Mx, My, Mz,  c='b', label='Milstein ', s=1)
        ax.legend(loc='best')
        plt.show()
        
        
        print(np.mean(np.array(Mx) - np.array(MPx)))
        print(np.mean(np.array(My) - np.array(MPy)))
        print(np.mean(np.array(Mz) - np.array(MPz)))
        
        
        return  #output
 

L = Lotka(h = 0.005 , N = 200)
L.test()