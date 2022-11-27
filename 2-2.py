from math import pi, sin
import numpy as np
import matplotlib.pyplot as plt

L = 1
N = 100
dx = L/(N+1)
X = [i*dx for i in range(N)]
def V(x):
    #return 0
    #return 700*(0.5-abs(x-0.5))
    return 800*sin(2*pi*x)

T = (np.diag(np.full(N,-2))+np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1))/(dx**2)
V_matrix = np.diag([V(x) for x in X])
T = T-V_matrix

print(T)
val, vectorsRow = np.linalg.eig(T)


vectors = [[] for _ in range(len(T))]
for i in range(len(T)):
    for v in vectorsRow:
        vectors[i].append(v[i])

val = list(val)

vectors = list(vectors)
sorted_vectors = [x for _, x in sorted(zip(val, vectors))]
sorted_val = sorted(val)


alpha = 0
beta = 0


dummy = ["first", "second", "thrid"]
plt.style.use('seaborn-poster')
plt.figure(figsize = (12, 8))


for i in range(3):
    eigenvector = sorted_vectors[-1-i]
    eigenvector = list(eigenvector)
    eigenvector.insert(alpha,0)
    eigenvector.append(beta)
    wave_func = [x for x in eigenvector]
    prob_dens = [sorted_val[-1-i]* x + abs(sorted_val[-1-i]) for x in wave_func]
    plt.plot([i*dx for i in range(len(prob_dens))], prob_dens, label=dummy[i])


plt.title('Wave function with energy levels for V(x) = 800*sin(2*pi*x)')
#plt.xlabel('t')
#plt.ylabel('f(t)')
plt.xlabel('x')
plt.ylabel('E_k psi(x) + E_k')
plt.grid()
plt.legend(loc='lower right')
plt.show()