from twopBVP import twopBVP, twopBVP_using_toeplitz
import math
import matplotlib.pyplot as plt
from numpy.linalg import norm
import numpy as np
different_N = [3**i for i in range(1,10)]
#N = 100
L = 10
E = 1.9*10**11
Qx = -50000
#dx = L/(N+1)

#X = [dx*i for i in range(1,N+1)] 
y0 = 0
yN_1 = L
α = 0
β = 0

Dxs = np.zeros(400-10)
rms = np.zeros(400-10)
def findGlobalError2():
    for N in range(10, 400):
        dx = L/(N+1)
        X = [dx*i for i in range(1,N+1)] 
        fvec = [Qx*dx**2 for x in X]
        x = np.linspace(0+dx,L-dx,N)
        y = twopBVP_using_toeplitz(fvec, 1, math.e**4, L, N)
        x = np.linspace(0, L, N+2)
        yReal = np.exp(np.square(x))
        err = np.zeros(len(y))
        for i in range(len(y)):
            err[i]=(abs(yReal[i] - y[i])**2*dx)
        rms[N-10] = math.sqrt(np.sum(err))
    plt.grid()
    plt.loglog(Dxs, rms)
    plt.show()

def find_global_error():
    globalError = []
    for N in different_N:
        dx = L/(N+1)
        Y = find_function(N)
        X = [dx*i for i in range(1,N+1)] 
        EXACT = [-math.sin(math.pi*x)/(math.pi**2) for x in X] # sin(pi*x)=y''
        #globalError.append(dx**0.5*norm(Y-EXACT))
        #globalError.append(dx**0.5*norm([Y[i-EXACT[i]] for i in range(len(Y))]))

        globalError.append(abs(Y[int(N/2)]-EXACT[int(N/2)]))
        #globalError.append(abs(Y[-3]-EXACT[-3]))

   
    plt.style.use('seaborn-poster')
    plt.figure(figsize = (12, 8))

    print(globalError)
    print([L/(i+1) for i in different_N])
    #plt.loglog([L/(i+1) for i in different_N], globalError, 'bo--', label='GlobalError')
    #plt.loglog([L/(i+1) for i in different_N], [0.06*(L/(i+1))**2 for i in different_N], 'g', label='x^2')
    plt.plot(X, Y, 'bo--', label='Approximate')
    plt.plot(X, EXACT, 'g', label='Exact')
    #plt.title('Approximate and Exact Solution for sin(pi*x)==y''')
    plt.title('Global error for sin(pi*x) for Beam equation')
    #plt.xlabel('t')
    #plt.ylabel('f(t)')
    plt.xlabel('N')
    plt.ylabel('GlobalError(N)')
    plt.grid()
    plt.legend(loc='lower right')
    plt.show()


def find_function(N):
    dx = L/(N+1)
    X = [dx*i for i in range(1,N+1)] 
    fvec = [Qx*dx**2 for x in X]
    fvec[0]-=α
    fvec[-1]-=β
    M = twopBVP_using_toeplitz(fvec, α, β , L, N)

    fvec = [M[i]/(E*f(dx*i))*dx**2 for i in range(N)]

    Y = twopBVP_using_toeplitz(fvec, 0, 0 , L, N)
    return Y
    # Y = M
    Y = list(Y)
    Y.insert(0,α)
    Y.append( β)

    X = list(X)
    X.insert(0,0)   
    X.append(L)




    plt.style.use('seaborn-poster')
    plt.figure(figsize = (12, 8))

    print([L/(i+1) for i in different_N])

    plt.plot(X, Y, 'bo--', label='Approximate')
    #plt.plot(X, EXACT, 'g', label='Exact')
    plt.title('Approximate and Exact Solution for Beam equation')

    plt.xlabel('t')
    plt.ylabel('f(t)')
    plt.grid()
    plt.legend(loc='lower right')
    plt.show()


def f(x):
    return 10**(-3)*(3-2*math.cos(math.pi*x/L)**12)
    #return math.exp(-x**2)

if __name__ == "__main__":
    #find_function()
    find_global_error()
