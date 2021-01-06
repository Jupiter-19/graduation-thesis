'''
Example 4.3
dX = sin(x) dt + sqrt(2) dB
dX = -x /sqrt(1+x**4) dt + sqrt(2) dB
可变参数：X0,T
精度参数：h,格式
'''
import numpy as np
import matplotlib.pyplot as plt

class Stable():
    def __init__(self,X0,h,T):
        self.X0 = X0
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
        X = [self.X0 for i in range(self.N+1)]

        for i in range(self.N):
            # Euler 格式 & Milstein 格式
            X[i+1] = X[i]  + np.sin(X[i]) * self.h + self.sqrt_2 * self.Brown[i]

        if myplot:
            plt.figure(figsize=(20,10))
            plt.plot(self.t,X)
            plt.show()
        
        return X[-1]
    
    def measure(self):
        result = []
        for idx in range(50000):
            if idx % 2000 == 0: print(idx//500)
            result.append(self.sample())

        plotRange = 100
        
        result = [x for x in result if x > -plotRange and x < plotRange]
        
        print(len(result))
        
        plt.figure(figsize=(18,9))
        plt.hist(result , bins = 250)
        plt.show()


T = 500
m = Stable(X0 = 0 , h = T/2500 , T = T)
#m.sample(myplot=1)
m.measure()




