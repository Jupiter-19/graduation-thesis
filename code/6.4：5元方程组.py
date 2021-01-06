
from scipy.optimize import root,fsolve
import numpy as np


a,b,c = 1,2,3
I1,I2 = 4,1

def f1(x):
    return np.array([x[0]-x[3]-x[4]*x[1]*x[2]-a,
                     x[1]-x[3]-x[4]*x[0]*x[2]-b,
                     x[2]-x[3]-x[4]*x[0]*x[1]-c,
                     x[0]+x[1]+x[2]-I1,
                     x[0]*x[1]*x[2]-I2])
 

sol_root = root(f1,[a,b,c,0,0])
print(sol_root.x)