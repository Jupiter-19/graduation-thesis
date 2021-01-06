'''
d X = [tan(-pi/2*x) + sign(x)] dt + |1-|x||^alpha X dw
方程参数: X0,b,sigma
时间参数：T
格式参数：h
'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class Method():
    def __init__(self):
        self.num = int(5e5)
        self.X0 =  [0.5] * (self.num//2) + [-0.5] * (self.num//2)
        self.alpha = 0.5
        
        
        self.h = 1e-3
        self.T = 2
        
        self.N = int(self.T / self.h)
        self.sqrt_h = np.sqrt(self.h)
        self.t = [self.h * i for i in range(self.N+1)]
            
        self.M1 = []
        self.M2 = []
    def sample(self , X):
        Brown = np.random.normal(0 , 1, self.num) * self.sqrt_h
        newX = []
        for idx,x in enumerate(X):
            a = np.tan(-np.pi/2*x) + np.sign(x)
            b = (1-np.abs(x)) ** self.alpha
            newx = x + a*self.h + b*Brown[idx]
            if newx >= 1:
                newX.append((x+1)/2)
            elif newx <= -1:
                newX.append((x-1)/2)
            else:
                newX.append(newx)
        
        return newX
    def test(self):
        X = self.X0
        for idx in range(self.N+1):
            newX = self.sample(X)
            tmp1 = np.sum(np.abs( np.sort(X)-np.sort(newX) )) / self.num
            tmp2 = np.sum( (np.sort(X)-np.sort(newX))**2 ) / self.num
            print(idx,tmp1,tmp2)
            self.M1.append(tmp1)
            self.M2.append(tmp2)
            X = newX
            if idx % 20 == 0:
                plt.figure(figsize=(16,8))
                plt.hist(newX,bins = 100)
                plt.show()
                plt.savefig('../images/6.2/measure1/' + str(idx) + '.png')
                out = pd.DataFrame({'M1':self.M1 , 'M2':self.M2})
                out.to_csv('./6.2.result2.csv')
        return X
# In[]
m = Method()
X = m.test()





