from qutip import *

# (1) K-space basis states, returns a list of all the basis states
def Kspace(N):
    L = []
    for i in range(N):
        L.append(qutip.basis(N,i))
    return L