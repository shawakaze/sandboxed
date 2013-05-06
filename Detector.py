from qutip import *
from Kspace import *
from Operators import *

def CPTPmap(A,B,p):
    a = A*p*A.dag() + B*p*B.dag()
    return a

def AntiCommutator(A,B):
    return A*B+B*A

def Lindblad(A,B,p):
    a = A*p*A.dag()+B*p*B.dag() - 0.5*AntiCommutator(A*A.dag()+B*B.dag(),p)
    return a

def Detection(A,B,p):
    Zero_ = tensor(Z(2),Z(2),Z(2),Z(9))
    g = Zero_
    for j in range(9):
        g = g + Lindblad(A[j],B[j],p)
    return g


