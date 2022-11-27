from scipy import sparse
import numpy as np
from scipy.linalg import solve_toeplitz

def twopBVP(fvec : np.array, α, β, L, N):
    h = L/(N+1)
    fvec[0] -= α
    fvec[N-1] -= β
    T = [[0 for i in range (N)] for i in range(N)]

    for i in range(N):
        T[i][i] = -2
        if i>0:
            T[i][i-1] = 1
        if i<N-1:
            T[i][i+1] = 1
    Tinv : np.array = sparse.linalg.inv(sparse.csc_matrix(T))
    return Tinv.dot(fvec)



def twopBVP_using_toeplitz(fvec : np.array, α, β, L, N):
    c = np.array([-2, 1] + [0 for i in range(N-2)])    # First column of T # First row of T
    b = fvec
    Y = solve_toeplitz((c, c), b)
    return Y


