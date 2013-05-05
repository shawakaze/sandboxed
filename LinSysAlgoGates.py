from qutip import *

import scipy as sp
import cmath

from Operators import *

"""
    (1): Ax=b
    (2): First we exponentiate the A matrix.. Luckily qutip has an exponentiation fucntion.
    (3): 
"""
# We define the global parameters here
i = cmath.sqrt(-1)
I = qeye(2)
pi = sp.pi

# the input matrix

A = 2*pi*i*matrix([[1,0],[0,1]])

print U

# (1) K-space basis states, returns a list of all the basis states
def Kspace(N):
    L = []
    for i in range(N):
        L.append(qutip.basis(N,i))
    return L