'''
d X = bX dt + sigma X dw
方程参数: X0,b,sigma
时间参数：T
格式参数：h
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class Method():
    def __init__(self,X0,b,sigma,h,T):
        self.X0 = X0
        self.b = b
        self.sigma = sigma
        self.h = h
        self.T = T
        self.set_para()
        print(h)
        
    def set_para(self):
        self.N = int(self.T / self.h)
        self.sqrt_h = np.sqrt(self.h)
        self.sqrt_3 = np.sqrt(3)
        self.t = [self.h * i for i in range(self.N+1)]
        
    def sample(self , myplot = 0):
        Brown = np.random.normal(0 , 1, self.N) * self.sqrt_h
        b , sigma , h = self.b , self.sigma ,self.h
        X = [self.X0 for i in range(self.N+1)]
        Y = [self.X0 for i in range(self.N+1)]
        Z = [self.X0 for i in range(self.N+1)]
        M = [self.X0 for i in range(self.N+1)]

        for i in range(self.N):
        	# Euler
            X[i+1] = X[i] + b*X[i]*h + sigma*X[i]*Brown[i]
            # Milstein
            Y[i+1] = Y[i] + b*Y[i]*h + sigma*Y[i]*Brown[i] + (Brown[i]**2-h)/2*Y[i]*sigma**2
            # Taylor
            Z[i+1] = Z[i] + b*Z[i]*h + sigma*Z[i]*Brown[i]
            Z[i+1] += (Brown[i]**2 - h)/2*Z[i]*sigma**2
            Z[i+1] += h**2/2 * b**2 * Z[i]
            Z[i+1] += Brown[i] * b * sigma * Z[i] * h
            Z[i+1] += (Brown[i]**3-3*h*Brown[i])/6 * sigma**3*Z[i]
            # midpoint
            tmp = b*h/2 - sigma**2*h/4 + sigma/2*Brown[i]
            M[i+1] = M[i]*(1+tmp)/(1-tmp)
        
        if myplot:
            plt.figure(figsize=(20,10))
            plt.plot(self.t,X)
        
        real = self.X0 * np.exp((b-sigma**2/2)*self.T + sigma*np.sum(Brown))
        
        sum_ = sum([X[-1] , Y[-1] , Z[-1] , M[-1] , real])
        if np.abs(sum_) > 1e6:
            print([X[-1] , Y[-1] , Z[-1] , M[-1] , real])
            print('error')
            
        return [X[-1] , Y[-1] , Z[-1] , M[-1] , real]
    
    def test(self):
        # EX = self.X0 * np.exp(self.b * self.T)
        # DX = (self.X0**2 * np.exp(2*self.b*self.T) * (np.exp(self.sigma**2*self.T)-1))
    	# print(EX,DX)
        
        euler = []
        milstein = []
        taylor = []
        midpoint = []
        for i in range(10000):
            if i %2000 == 0:
                print(i)
            calc1,calc2,calc3,calc4,real = self.sample() 
            euler.append(np.abs(calc1 - real))
            milstein.append(np.abs(calc2 - real))
            taylor.append(np.abs(calc3 - real))
            midpoint.append(np.abs(calc4 - real))
        m_mean = [np.mean(euler),np.mean(milstein),np.mean(taylor),np.mean(midpoint)]
        m_std  = [np.std(euler),np.std(milstein),np.std(taylor),np.std(midpoint)]

        print(m_mean)
        print(m_std)
        
        return m_mean , m_std

       
# In[]
Euler1 = []
Milstein1 = []
Taylor1 = []
Midpoint1 = []

Euler2 = []
Milstein2 = []
Taylor2 = []
Midpoint2 = []



hlim = [1e-5,2e-5,5e-5,1e-4,2e-4,5e-4,1e-3]
for idx in range(7):
    m = Method(X0 = 2 , b = 2, sigma = 5 , h = hlim[idx], T = 0.1)
    M,S = m.test()
    
    Euler1.append(M[0])
    Milstein1.append(M[1])
    Taylor1.append(M[2])
    Midpoint1.append(M[3])
    
    Euler2.append(S[0])
    Milstein2.append(S[1])
    Taylor2.append(S[2])
    Midpoint2.append(S[3])

Dict = {'Euler1':Euler1 , 'Milstein1':Milstein1 , 'Taylor1':Taylor1 , 'Midpoint1':Midpoint1,
        'Euler2':Euler2 , 'Milstein2':Milstein2 , 'Taylor2':Taylor2 , 'Midpoint2':Midpoint2,}
out = pd.DataFrame(Dict)
out.to_csv('6.1.result.csv')

# In[]
hlim = [1e-5,2e-5,5e-5,1e-4,2e-4,5e-4,1e-3]
data = pd.read_csv('6.1.result.csv')

Euler1 = list(data['Euler1'])
Euler2 = list(data['Euler2'])

Milstein1 = list(data['Milstein1'])
Milstein2 = list(data['Milstein2'])

Taylor1 = list(data['Taylor1'])
Taylor2 = list(data['Taylor2'])

Midpoint1 = list(data['Midpoint1'])
Midpoint2 = list(data['Midpoint2'])

t = [np.log(np.abs(hlim[idx])) for idx in range(7)]
p1= [np.log(np.abs(Euler1[idx])) for idx in range(7)]
p2= [np.log(np.abs(Milstein1[idx])) for idx in range(7)]
p3= [np.log(np.abs(Taylor1[idx])) for idx in range(7)]
p4= [np.log(np.abs(Midpoint1[idx])) for idx in range(7)]

plt.figure(figsize=(16,8))
plt.xlabel('log(h)')
plt.ylabel('log(E(error))')
plt.plot(t,p1,label='Euler',    color='r',marker='o')
plt.plot(t,p2,label='Milstein', color='g',marker='v')
plt.plot(t,p3,label='Taylor',   color='y',marker='*')
plt.plot(t,p4,label='Midpoint', color='b',marker='s')
plt.legend() #显示图例
plt.show()

    
t = [np.log(np.abs(hlim[idx])) for idx in range(7)]
p1= [np.log(np.abs(Euler2[idx])**2) for idx in range(7)]
p2= [np.log(np.abs(Milstein2[idx])**2) for idx in range(7)]
p3= [np.log(np.abs(Taylor2[idx])**2) for idx in range(7)]
p4= [np.log(np.abs(Midpoint2[idx])**2) for idx in range(7)]

plt.figure(figsize=(16,8))
plt.xlabel('log(h)')
plt.ylabel('log(E(error**2))')
plt.plot(t,p1,label='Euler',    color='r',marker='o')
plt.plot(t,p2,label='Milstein', color='g',marker='v')
plt.plot(t,p3,label='Taylor',   color='y',marker='*')
plt.plot(t,p4,label='Midpoint', color='b',marker='s')
plt.legend() #显示图例
plt.show()

    
    
    


