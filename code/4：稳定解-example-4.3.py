'''
方程：dX = bx dt + sqrt(2)*(1-x**2) dB
固定参数：b
可变参数：T
精度参数：h,格式
'''
import numpy as np
import matplotlib.pyplot as plt

class Stable():
    def __init__(self,b,h,T):
        self.b = b
        self.h = h
        self.T = T
        self.set_para()

    def set_para(self):
        self.N = int(self.T / self.h)
        self.sqrt_h = np.sqrt(self.h)
        self.sqrt_2 = np.sqrt(2)
        self.t = [self.h * i for i in range(self.N+1)]

    def sample(self,myplot = 0):
        self.Brown = np.random.normal(0 , 1, self.N) * self.sqrt_h
        X = [0 for i in range(self.N+1)]
        
        tmp = np.random.normal(0,0.2)
        while np.abs(tmp) >= 1:
            tmp = np.random.normal(0,0.2)
        X[0] = tmp        
        

        for i in range(self.N):
            # Milstein 格式
            X[i+1] = X[i] + self.b * X[i] * self.h + self.sqrt_2 * (1-X[i]**2) * self.Brown[i]
            X[i+1]+= (self.Brown[i]**2 - self.h) * 2 * X[i] * (X[i]**2-1)
            
            if X[i+1] > 1: X[i+1] = (X[i]+1) / 2
            if X[i+1] <-1: X[i+1] = (X[i]-1) / 2
        
        if myplot:
            plt.figure(figsize=(20,10))
            plt.plot(self.t,X)
        
        return X[-1]
    
    def measure(self):
        result = []
        for idx in range(500000):
            if idx % 10000 == 0: print(idx//5000+2)
            result.append(self.sample())

        plt.figure(figsize=(18,9))
        plt.xlim(-1,1)
        
        plt.hist(result , bins = 250)
        plt.show()


T = 20.0
m = Stable(b = -2 , h = T/1600 , T = T)
#m.sample(myplot=1)
m.measure()








