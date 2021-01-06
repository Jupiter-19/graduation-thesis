import numpy as np
import matplotlib.pyplot as plt

class Kubo():
    def __init__(self,h):
        self.X0, self.Y0 = 1 , 0
        self.b, self.sigma = 1 , 1
        self.h = h
        self.T = 20
        
        
        self.num = 5000
        self.N = int(self.T / self.h)
        
        self.sqrt_h = np.sqrt(self.h)
        self.t = [self.h * i for i in range(self.N+1)]


    def scheme(self, Ex,Ey,Mx,My,Tx,Ty,EPx,EPy,MPx,MPy ):
        b,sigma,h = self.b , self.sigma, self.h
        Brown = np.random.normal(0 , 1, self.num) * self.sqrt_h
        
        ex = [0 for i in range(self.num)]
        ey = [0 for i in range(self.num)]
        mx = [0 for i in range(self.num)]
        my = [0 for i in range(self.num)]
        tx = [0 for i in range(self.num)]
        ty = [0 for i in range(self.num)]
        
        epx= [0 for i in range(self.num)]
        epy= [0 for i in range(self.num)]
        mpx= [0 for i in range(self.num)]
        mpy= [0 for i in range(self.num)]
        
        
    
        
        for i in range(self.num):
            # Euler格式 Euler-Maruyama method
            ex[i] = Ex[i] - b*Ey[i]*h - sigma*Ey[i]*Brown[i] - sigma**2*Ex[i]/2*h
            ey[i] = Ey[i] + b*Ex[i]*h + sigma*Ex[i]*Brown[i] - sigma**2*Ey[i]/2*h
            # Milstein格式
            mx[i] = Mx[i] - b*My[i]*h - sigma*My[i]*Brown[i] - sigma**2*Mx[i]/2*Brown[i]**2
            my[i] = My[i] + b*Mx[i]*h + sigma*Mx[i]*Brown[i] - sigma**2*My[i]/2*Brown[i]**2
            # Taylor格式
            tx[i] = Tx[i] - b*Ty[i]*h - sigma*Ty[i]*Brown[i] - sigma**2*Tx[i]/2*Brown[i]**2
            tx[i] +=(Brown[i]**3-3*h*Brown[i])/6 * sigma**3*Ty[i]
            tx[i] +=(sigma**4/4*Tx[i]-sigma**2*b*Ty[i]-b**2*Tx[i])/2*h*h
            tx[i] -= Brown[i]*h*(sigma**3/2*Ty[i]+b*sigma*Tx[i])
            ty[i] = Ty[i] + b*Tx[i]*h + sigma*Tx[i]*Brown[i] - sigma**2*Ty[i]/2*Brown[i]**2
            ty[i] -=(Brown[i]**3-3*h*Brown[i])/6 * sigma**3*Tx[i]
            ty[i] +=(sigma**4/4*Ty[i]+sigma**2*b*Tx[i]-b**2*Ty[i])/2*h*h
            ty[i] += Brown[i]*h*(sigma**3/2*Tx[i]-b*sigma*Ty[i])
            # EulerP           
            epx[i] = EPx[i] - b*EPy[i]*h - sigma*EPy[i]*Brown[i] - sigma**2*EPx[i]/2*h
            epy[i] = EPy[i] + b*EPx[i]*h + sigma*EPx[i]*Brown[i] - sigma**2*EPy[i]/2*h            
            tmp = np.sqrt(epx[i]**2 + epy[i]**2)
            epx[i] = epx[i] / tmp
            epy[i] = epy[i] / tmp
            # MilsteinP
            mpx[i] = MPx[i] - b*MPy[i]*h - sigma*MPy[i]*Brown[i] - sigma**2*MPx[i]/2*Brown[i]**2
            mpy[i] = MPy[i] + b*MPx[i]*h + sigma*MPx[i]*Brown[i] - sigma**2*MPy[i]/2*Brown[i]**2
            tmp = np.sqrt(mpx[i]**2 + mpy[i]**2)
            mpx[i] = mpx[i] / tmp
            mpy[i] = mpy[i] / tmp
            
        return [ex,ey,mx,my,tx,ty,epx,epy,mpx,mpy]
            
    
    

    def test(self):
        Ex, Ey = [1]*self.num , [0]*self.num # Euler 格式
        Mx, My = [1]*self.num , [0]*self.num # Milstein 格式
        Tx, Ty = [1]*self.num , [0]*self.num # Taylor 格式
        EPx,EPy= [1]*self.num , [0]*self.num # Euler 格式
        MPx,MPy= [1]*self.num , [0]*self.num # Milstein 格式
        
        for idx in range(self.N):
            print(round(self.h*idx,4))
            output = self.scheme(Ex,Ey,Mx,My,Tx,Ty,EPx,EPy,MPx,MPy)
            Ex,Ey,Mx,My,Tx,Ty,EPx,EPy,MPx,MPy = output
        
        return  output
 

k = Kubo(h = 0.01)
Ex,Ey,Mx,My,Tx,Ty,EPx,EPy,MPx,MP = k.test()

# In[]
plt.figure(figsize=(12,12))

ex,ey,mx,my,epx,epy = [],[],[],[],[],[]
for i in range(5000):
    if np.abs(Ex[i]) < 1.5 and np.abs(Ey[i]) < 1.5:
        ex.append(Ex[i])
        ey.append(Ey[i])
    if np.abs(Mx[i]) < 1.5 and np.abs(My[i]) < 1.5:
        mx.append(Mx[i])
        my.append(My[i])

plt.scatter(ex,ey, label='Euler',    color='y', s=2,marker='v')
plt.scatter(mx,my, label='Milstein', color='r', s=2,marker='*')
plt.scatter(EPx,EPy, label='EulerP', color='b', s=2,marker='s')
plt.legend()
plt.xlim=[-1,1]
plt.show()


np.sum(np.array(Ex) - np.array(Tx))