import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import root,fsolve
from mpl_toolkits.mplot3d import Axes3D 

class Lotka():
    def __init__(self,h,N,X,Y,Z):
        self.h = h
        self.N = N
        self.I1 = 4
        self.I2 = 2
        self.X = X
        self.Y = Y
        self.Z = Z
        
        self.num = len(self.X) # 样本轨道数

    def loss(self,x):
        return np.array([x[0]-x[3]-x[4]*x[1]*x[2]-self.a,
                         x[1]-x[3]-x[4]*x[0]*x[2]-self.b,
                         x[2]-x[3]-x[4]*x[0]*x[1]-self.c,
                         x[0]+x[1]+x[2]-self.I1,
                         x[0]*x[1]*x[2]-self.I2])
    
    
    def scheme(self, MPx,MPy,MPz):
        h = self.h
        Brown = np.random.normal(0 , 1, self.num) * np.sqrt(h)
        
        mpx= [0 for i in range(self.num)]
        mpy= [0 for i in range(self.num)]
        mpz= [0 for i in range(self.num)]
        
        
    
        for i in range(self.num):
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

        return [mpx,mpy,mpz]
            
    
    

    def test(self):
        MPx,MPy,MPz = self.X , self.Y, self.Z # MilsteinP 格式
        
        for idx in range(self.N):
            print(round(self.h*idx,4))
            MPx,MPy,MPz = self.scheme(MPx,MPy,MPz)

        plt.rcParams['font.sans-serif']=['SimHei'] # 显示中文
        plt.rcParams['axes.unicode_minus']=False #用来正常显示负号
        plt.figure(figsize=(16,8))
        plt.title('X的稳定解')
        plt.hist(MPx , bins = 100)
        plt.show()
        
        plt.figure(figsize=(16,8))
        plt.title('Y的稳定解')
        plt.hist(MPy , bins = 100)
        plt.show()
        
        plt.figure(figsize=(16,8))
        plt.title('Z的稳定解')
        plt.hist(MPz , bins = 100)
        plt.show()
        
        
        out = pd.DataFrame({'X':MPx , 'Y':MPy , 'Z':MPz})
        out.to_csv('6.4.stable.csv')
        
        I1 = [MPx[i]+MPy[i]+MPz[i] for i in range(self.num)]
        I2 = [MPx[i]*MPy[i]*MPz[i] for i in range(self.num)]
        
        print(np.mean(I1) , np.std(I1))
        print(np.mean(I2) , np.std(I2))
        
        print(self.num)
        
        return  #output
 
data = pd.read_csv('6.4.stable.csv')
X,Y,Z = list(data['X']), list(data['Y']) , list(data['Z'])

L = Lotka(h = 0.001 , N = 5 , X = X, Y = Y, Z = Z)
L.test()










