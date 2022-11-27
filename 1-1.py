from twopBVP import twopBVP, twopBVP_using_toeplitz
import math
import matplotlib.pyplot as plt
from numpy.linalg import norm

different_N = [3**i for i in range(1,10)]

def f(x):
    return math.sin(math.pi*x) # sin(pi*x)
    #return math.exp(-x**2)

if __name__ == "__main__":
    globalError = []
    for N in different_N:
        L = 1
        dx = L/(N+1)
        X = [dx*i for i in range(1,N+1)] 
        y0 = 0
        yN_1 = 1
        α = 0
        β = 0
        const = dx**2
        fvec = [f(x)*const for x  in X]
        fvec[0]-=α
        fvec[-1]-=β


        Y = twopBVP_using_toeplitz(fvec, α, β , L, N)
        EXACT = [-math.sin(math.pi*x)/(math.pi**2) for x in X] # sin(pi*x)=y''
        globalError.append(norm(Y-EXACT))
        Y = list(Y)
        Y.insert(0,α)
        Y.append( β)

        X = list(X)
        X.insert(0,0)   
        X.append(L)
        EXACT = [-math.sin(math.pi*x)/(math.pi**2) for x in X] # sin(pi*x)=y''
      
        #globalError.append(abs(Y[int(N/2)]-EXACT[int(N/2)]))
        #globalError.append(abs(Y[-3]-EXACT[-3]))


    plt.style.use('seaborn-poster')
    plt.figure(figsize = (12, 8))

    #print(globalError)
    #print([L/(i+1) for i in different_N])
    plt.loglog([ L/(i+1) for i in different_N], globalError, 'bo--', label='GlobalError')
    plt.loglog([ L/(i+1) for i in different_N], [0.06*(L/(i+1))**2 for i in different_N], 'g', label='x^2')
    #plt.plot(X, Y, 'bo--', label='Approximate')
    #plt.plot(X, EXACT, 'g', label='Exact')
    #plt.title('Approximate and Exact Solution for sin(pi*x)==y''')
    plt.title('Global error for sin(pi*x) for different N')
    #plt.xlabel('t')
    #plt.ylabel('f(t)')
    plt.xlabel('N')
    plt.ylabel('GlobalError(N)')
    plt.grid()
    plt.legend(loc='lower right')
    plt.show()



