import numpy as np
import matplotlib.pyplot as plt
N = 100
L = 1
dx = L/(N+1)
const = 1/dx**2
#c = np.array([-2, 1] + [0 for i in range(N-2)])    # First column of T # First row of T

T = [[0 for i in range (N)] for i in range(N)]


for i in range(N):
    T[i][i] = -2*const
    if i>0:
        T[i][i-1] = 1*const
    if i<N-1:
        T[i][i+1] = 1*const

T[-1][-2] = 2*const

#for line in T:
#    print(line)


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
    yn = 2*beta+eigenvector[-2]
    eigenvector.append(yn)
    eigenvector = [x*yn for x in eigenvector]
    plt.plot([i*dx for i in range(len(eigenvector))], eigenvector, label=dummy[i])
    

plt.title('Eigenfunctions of d^2/dx^2')
#plt.xlabel('t')
#plt.ylabel('f(t)')
plt.xlabel('x')
plt.ylabel('eigenvector values')
plt.grid()
plt.legend(loc='lower right')
plt.show()