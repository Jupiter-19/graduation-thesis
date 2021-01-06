import numpy as np
import matplotlib.pyplot as plt

N = 200

t = [i/N*2 for i in range(-N//2,N//2+1)]

y = []
for i in range(N):
    gap = (np.arcsin(t[i+1]) - np.arcsin(t[i]) ) / np.pi
    prob = round(gap * 2e5)
    y += list(np.random.uniform(t[i],t[i+1],int(prob)))
    
plt.figure(figsize=(16,8))
plt.hist(y,bins = N//2)