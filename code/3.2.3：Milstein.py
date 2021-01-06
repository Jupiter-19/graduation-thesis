'''
d X = bX dt + sigma X dw
方程参数: X0,b,sigma
时间参数：T
格式参数：h
'''
import numpy as np
import matplotlib.pyplot as plt

class Method():
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
        self.sqrt_3 = np.sqrt(3)
        self.t = [self.h * i for i in range(self.N+1)]
        
    def sample(self , myplot = 0):
        self.Brown = np.random.normal(0 , 1, self.N) * self.sqrt_h
        #self.Brown2= np.random.normal(0 , 1, self.N) * self.sqrt_h
        X = [self.X0 for i in range(self.N+1)]
        Y = [self.X0 for i in range(self.N+1)]
        Z = [self.X0 for i in range(self.N+1)]

        for i in range(self.N):
        	# Euler
            #X[i+1] = X[i] + self.b * X[i] * self.h + self.sigma * X[i] * self.Brown[i]
            # Euler2
            X[i+1] = X[i] + self.b * X[i] * self.h + self.sigma * X[i] * self.Brown[i]
            X[i+1] =self.sigma**2/2 * X[i] * self.h
            # Milstein
            Y[i+1] = Y[i] + self.b * Y[i] * self.h + self.sigma * Y[i] * self.Brown[i]
            Y[i+1] += (self.Brown[i]**2 - self.h)/2 * Y[i] * self.sigma**2
            # Taylor
            Z[i+1] = Z[i] + self.b * Z[i] * self.h + self.sigma * Z[i] * self.Brown[i]
            Z[i+1] += (self.Brown[i]**2 - self.h)/2 * Z[i] * self.sigma**2
            Z[i+1] += self.h**2/2 * self.b**2 * Z[i]
            Z[i+1] += self.Brown[i] * self.b * self.sigma * Z[i] * self.h
            Z[i+1] += (self.Brown[i]**3-3*self.h*self.Brown[i])/6 * self.sigma**3*Z[i]
        
        if myplot:
            plt.figure(figsize=(20,10))
            plt.plot(self.t,X)
        
        real = self.X0 * np.exp((self.b-self.sigma**2/2)*self.T+self.sigma*np.sum(self.Brown))
        return (X[-1] , Y[-1] , Z[-1] , real)
    
    def test(self):
        # EX = self.X0 * np.exp(self.b * self.T)
        # DX = (self.X0**2 * np.exp(2*self.b*self.T) * (np.exp(self.sigma**2*self.T)-1))
    	# print(EX,DX)
        
        euler = []
        milstein = []
        taylor = []
        for i in range(10000):
            if i %2000 == 0:
                print(i)
            calc1,calc2,calc3,real = self.sample()
            euler.append(calc1 - real)
            milstein.append(calc2 - real)
            taylor.append(calc3 - real)
        m_mean = [np.mean(euler),np.mean(milstein),np.mean(taylor)]
        m_std  = [np.std(euler),np.std(milstein),np.std(taylor)]

        print(m_mean)
        print(m_std)

       
m = Method(X0 = 2 , b = 2, sigma = 5 , h = 1e-2, T = 0.5)
#m.sample(myplot = 1)
m.test()








