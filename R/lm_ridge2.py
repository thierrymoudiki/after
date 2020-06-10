import numpy as np 
from scipy.linalg import lstsq
from time import time 
import matplotlib.pyplot as plt


n = 100
p = 10
lam = 0.1
lams = 10**np.linspace(start=-2, stop=3, num=25)
betas = []


y = np.concatenate((np.random.rand(n), np.repeat(0, p)))

for lam in lams:
    X = np.vstack((np.random.rand(n, p), np.sqrt(lam)*np.eye(p)))
    betas.append(np.linalg.lstsq(X, y, rcond=None)[0])


plt.plot(np.asarray(betas))
plt.show()

print(lams)

print(np.log(lams))
