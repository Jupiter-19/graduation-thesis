import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class Kubo():
    def __init__(self,h,X,Y,N):
        self.b, self.sigma = 1 , 1
        self.h = h
        self.N = N
        self.X = X
        self.Y = Y
        self.num = len(X)
        self.sqrt_h = np.sqrt(self.h)


    def scheme(self, Tx,Ty ):
        b,sigma,h = self.b , self.sigma, self.h
        Brown = np.random.normal(0 , 1, self.num) * self.sqrt_h
        tx = [0 for i in range(self.num)]
        ty = [0 for i in range(self.num)]
        
        for i in range(self.num):
            # Taylor格式
            tx[i] = Tx[i] - b*Ty[i]*h - sigma*Ty[i]*Brown[i] - sigma**2*Tx[i]/2*Brown[i]**2
            tx[i] +=(Brown[i]**3-3*h*Brown[i])/6 * sigma**3*Ty[i]
            tx[i] +=(sigma**4/4*Tx[i]-sigma**2*b*Ty[i]-b**2*Tx[i])/2*h*h
            tx[i] -= Brown[i]*h*(sigma**3/2*Ty[i]+b*sigma*Tx[i])
            ty[i] = Ty[i] + b*Tx[i]*h + sigma*Tx[i]*Brown[i] - sigma**2*Ty[i]/2*Brown[i]**2
            ty[i] -=(Brown[i]**3-3*h*Brown[i])/6 * sigma**3*Tx[i]
            ty[i] +=(sigma**4/4*Ty[i]+sigma**2*b*Tx[i]-b**2*Ty[i])/2*h*h
            ty[i] += Brown[i]*h*(sigma**3/2*Tx[i]-b*sigma*Ty[i])
            # TaylorP
            tmp1 = np.sqrt(Tx[i]**2 + Ty[i]**2) 
            tmp2 = np.sqrt(tx[i]**2 + ty[i]**2) 
            tx[i] = tx[i] / tmp2 * tmp1
            ty[i] = ty[i] / tmp2 * tmp1
            
        return tx , ty
            
    
    

    def test(self):
        Tx, Ty = self.X , self.Y # Taylor 格式
        
        for idx in range(self.N):
            print(round(self.h*idx,4))
            Tx,Ty = self.scheme(Tx,Ty)
        
        return Tx,Ty
 

# In[]

data = pd.read_csv('6.3.stable.csv')
X = list(data["X"])
Y = list(data["Y"])

#X = [(1+i)/40000 for i in range(39999)]
#Y = [(1+i)/40000 for i in range(39999)]

#X = X*2
#Y = Y*2

h = 0.001
N = 10

print(h,N)
k = Kubo(h, X, Y , N)
X,Y = k.test()
plt.rcParams['font.sans-serif']=['SimHei'] # 显示中文
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号
plt.figure(figsize=(16,8))
plt.title('X的稳定解')
plt.hist(X , bins = 200)
plt.show()

plt.figure(figsize=(16,8))
plt.title('Y的稳定解')
plt.hist(Y , bins = 200)
plt.show()

I = [X[i]**2 + Y[i]**2 for i in range(len(X))]

print(np.mean(I) , np.std(I) , len(I))

out = pd.DataFrame({"X":X,"Y":Y})
out.to_csv('6.3.stable.csv')





