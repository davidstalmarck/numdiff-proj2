import numpy as np
from math import pi
import matplotlib.pyplot as plt
L = 1
minN = 10
maxN=499

errors = [[],[],[]]
for n in range(minN, maxN):
    dx = L/(n+1)
    const = 1/dx**2
    #c = np.array([-2, 1] + [0 for i in range(N-2)])    # First column of T # First row of T

    #T = [[0 for i in range (n)] for i in range(n)]


    #for i in range(n):
    #    T[i][i] = -2*const
    #    if i>0:
    #        T[i][i-1] = 1*const
    #    if i<n-1:
    #        T[i][i+1] = 1*const

    #T[-1][-2] = 2*const
    T = np.diag(np.full(n,-2*const))+np.diag(np.ones(n-1)*const,1)+np.diag(np.ones(n-1)*const,-1)
    T[-1][-2] = 2*const

    val, vectors = np.linalg.eig(T)

    val = list(val)
    sorted_val = sorted(val)


    errors[0].append(abs(sorted_val[-1] + (pi/2)**2))
    errors[1].append(abs(sorted_val[-2] + (3*pi/2)**2))
    errors[2].append(abs(sorted_val[-3] + (5*pi/2)**2))
    
dummy = ["first", "second", "thrid"]
plt.style.use('seaborn-poster')
plt.figure(figsize = (12, 8))




for i in range(3):
    plt.loglog(list(range(minN,maxN)), errors[i], label=dummy[i])

plt.loglog(list(range(minN,maxN)),[10*i**(-1) for i in list(range(minN,maxN))], label="10*x^-1 in loglog")

plt.title('Eigenvalue global error')
#plt.xlabel('t')
#plt.ylabel('f(t)')
plt.xlabel('dx')
plt.ylabel('globalError(dx)')
plt.grid()
plt.legend(loc='lower right')
plt.show()

