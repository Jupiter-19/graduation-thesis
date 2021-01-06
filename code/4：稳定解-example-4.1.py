'''
Example 4.1
dX = bx dt + sqrt(2*sigma*(x**2+1)) dB
固定参数：b,sigma
可变参数：X0,T
精度参数：h,格式
'''
import numpy as np
import matplotlib.pyplot as plt

class Stable():
    def __init__(self,X0,b,sigma,h,T):
        self.X0 = X0
        self.b = b
        self.sigma = sigma
        self.h = h
        self.T = T
        self.set_para()

    def set_para(self):
        self.N = int(self.T / self.h)
        self.sqrt_h = np.sqrt(self.h)
        self.t = [self.h * i for i in range(self.N+1)]

    def sample(self,myplot = 0):
        self.Brown = np.random.normal(0 , 1, self.N) * self.sqrt_h
        X = [self.X0 for i in range(self.N+1)]
        Y = [self.X0 for i in range(self.N+1)]

        for i in range(self.N):
            # Euler 格式
            X[i+1] = X[i] + self.b * X[i] * self.h + np.sqrt(2*self.sigma*(X[i]**2+1)) * self.Brown[i]
            # Milstein 格式
            Y[i+1] = Y[i] + self.b * Y[i] * self.h + np.sqrt(2*self.sigma*(Y[i]**2+1)) * self.Brown[i]
            Y[i+1]+= (self.Brown[i]**2 - self.h) * Y[i] * self.sigma
        
        if myplot:
            plt.figure(figsize=(20,10))
            plt.plot(self.t,X)
            plt.show()
        
        return X[-1],Y[-1]
    
    def measure(self):
        result1,result2,error = [],[],[]
        for idx in range(50000):
            if idx % 500 == 0: print(idx)
            output1,output2 = self.sample()
            result1.append(output1)
            result2.append(output2)
            error.append(output1-output2)

        plotRange = 20
        
        result1 = [x for x in result1 if x > -plotRange and x < plotRange]
        result2 = [x for x in result2 if x > -plotRange and x < plotRange]
        
        # print(len(result1) , len(result2))
        # print(np.mean(error) , np.std(error))
        
        #plt.hist(result1 , bins = 100)
        #plt.show()
        plt.figure(figsize=(18,9))
        plt.xlim(-plotRange,plotRange)
        plt.hist(result2 , bins = 500)
        plt.show()

T = 100
m = Stable(X0 = 10 , b = 2, sigma = 5 , h = 0.1 , T = T)
#m.sample(myplot=1)
m.measure()
