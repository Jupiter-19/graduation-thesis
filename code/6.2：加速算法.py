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
            
        self.M1 = []
        self.M2 = []
        
    def sample(self , X, h):
        Brown = np.random.normal(0 , 1, self.num) * np.sqrt(h)
        newX = []
        for idx,x in enumerate(X):
            a = np.tan(-np.pi/2*x) + np.sign(x)
            b = (1-np.abs(x)) ** self.alpha
            newx = x + a*h + b*Brown[idx]
            if newx >= 1:
                newX.append((x+1)/2)
            elif newx <= -1:
                newX.append((x-1)/2)
            else:
                newX.append(newx)
        
        return newX
    
    def test(self):
        X = self.X0
        h = 0.1
        
        flag = 0
        for i in range(10):
            h /= 2
            M = int(np.power(2,i))
            for j in range(M):
                newX = self.sample(X,h)
                tmp1 = np.sum(np.abs( np.sort(X)-np.sort(newX) )) / self.num
                tmp2 = np.sum( (np.sort(X)-np.sort(newX))**2 ) / self.num
                print(flag,i,j,h,tmp1,tmp2)
                self.M1.append(tmp1)
                self.M2.append(tmp2)
                X = newX
                if flag % 20 == 0:
                    plt.figure(figsize=(16,8))
                    plt.hist(newX,bins = 100)
                    plt.savefig('C:/Users/Jupiter/Desktop/我的坚果云/研三上/.毕业设计/images/6.2/measure/' + str(i) 
                                + '.' + str(j) + '.png')
                    plt.show()
                    out = pd.DataFrame({'M1':self.M1 , 'M2':self.M2})
                    out.to_csv('./6.2.result2.csv')
                flag += 1
        return X
m = Method()
X = m.test()





