import numpy as np
import matplotlib.pyplot as plt

class Kubo():
    def __init__(self, init_value, paras, h, T):
        self.X0, self.Y0 = init_value
        self.b, self.sigma = paras 
        self.h = h
        self.T = T
        self.set_para()
        
    def set_para(self):
        self.N = int(self.T / self.h)
        self.sqrt_h = np.sqrt(self.h)
        self.t = [self.h * i for i in range(self.N+1)]
        self.sqrt_I = np.sqrt(self.X0**2 + self.Y0**2)

    def Euler(self,myplot = 0):
        Brown = np.random.normal(0 , 1, self.N) * self.sqrt_h
        X = [self.X0 for i in range(self.N+1)]
        Y = [self.Y0 for i in range(self.N+1)]
        
        for i in range(self.N):
            X[i+1] = X[i] - self.b * Y[i] * self.h - self.sigma * Y[i] * Brown[i]
            Y[i+1] = Y[i] + self.b * X[i] * self.h + self.sigma * X[i] * Brown[i]
            
        if myplot:
            plt.figure(figsize=(8,8))
            plt.plot(X,Y)
            plt.axis([-1,1,-1,1])
            plt.show()
        return X[-1],Y[-1]
    
    def EulerP(self,myplot = 0):
        Brown = np.random.normal(0 , 1, self.N) * self.sqrt_h
        X = [self.X0 for i in range(self.N+1)]
        Y = [self.Y0 for i in range(self.N+1)]
        
        for i in range(self.N):
            X[i+1] = X[i] - self.b * Y[i] * self.h - self.sigma * Y[i] * Brown[i]
            X[i+1]+= self.sigma**2 / 2 * Brown[i] ** 2 * X[i]
            Y[i+1] = Y[i] + self.b * X[i] * self.h + self.sigma * X[i] * Brown[i]
            Y[i+1]+= self.sigma**2 / 2 * Brown[i] ** 2 * Y[i]
            tmp = np.sqrt(X[i+1]**2 + Y[i+1]**2)
            X[i+1] = X[i+1] / tmp * self.sqrt_I
            Y[i+1] = Y[i+1] / tmp * self.sqrt_I
            
        if myplot:
            plt.figure(figsize=(8,8))
            plt.scatter(X,Y)
            plt.show()
        return X[-1],Y[-1]
            
    def stable(self):
        output1,output2 = [],[]
        for idx in range(int(5e4)):
            if idx % 2000 == 0: print(idx//500)
            x,y = self.EulerP()
            output1.append(x)
            output2.append(y)
        valueRange = 5
        output1 = [x for x in output1 if x < valueRange and x > - valueRange]
        output2 = [y for y in output2 if y < valueRange and y > - valueRange]
        
        plt.figure(figsize=(16,8))
        plt.hist(output1,bins = 250)
        plt.show()
        plt.figure(figsize=(16,8))
        plt.hist(output2,bins = 250)
        plt.show()
        
        plt.figure(figsize=(8,8))
        plt.axis([-1,1,-1,1])
        plt.scatter(output1,output2)
        plt.show()
        
        x = [i/10000 for i in range(-9900,9900)]
        y = [np.abs(t) / np.sqrt(1-t**2) for t in x]
        plt.plot(x,y)
    
k = Kubo(init_value=(1,0) , paras=(2,2) , h=1e-1, T=10)
#k.EulerP(myplot=1)
k.stable()


# Euler 方法下，无稳定解
# EulerP 方法下，有稳定解





