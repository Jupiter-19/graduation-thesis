'''
d X = bX dt + sigma X dw
方程参数: X0,b,sigma
时间参数：T
格式参数：h
'''
import numpy as np
import matplotlib.pyplot as plt

class Euler():
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
        
    def sample(self , myplot = 0):
        self.Brown = np.random.normal(0 , 1, self.N)
        self.Brown *= self.sqrt_h
        X = [self.X0 for i in range(self.N+1)]

        for i in range(self.N):
            X[i+1] = X[i] + self.b * X[i] * self.h + self.sigma * X[i] * self.Brown[i]
        
        if myplot:
            plt.figure(figsize=(20,10))
            plt.plot(self.t,X)
        
        real = self.X0 * np.exp((self.b-self.sigma**2/2)*self.T+self.sigma*np.sum(self.Brown))
        return (X[-1] , real)
    
    def test(self):
        EX = self.X0 * np.exp(self.b * self.T)
        DX = (self.X0**2 * np.exp(2*self.b*self.T) * (np.exp(self.sigma**2*self.T)-1))
    
        print(EX,DX)
        
        result = []
        error = []
        for i in range(10000):
            if i %2000 == 0:
                print(i)
            calc,real = self.sample()
            result.append(calc)
            error.append(real-calc)
        m_mean = np.mean(result) 
        m_std  = np.std(result) 
        print(m_mean , m_std ** 2)

        print(np.mean(error) , np.std(error)**2)

        
        
        result = [ x for x in result if x > m_mean-5*m_std and x < m_mean + 5*m_std]
        
        
        plt.hist(result,bins = 100)
        
m = Euler(X0 = 2 , b = 2, sigma = 5 , h = 1e-1 , T = 2)
#m.sample(myplot = 1)
m.test()








