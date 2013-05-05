from qutip import *
from UnitaryGates import *
from cmath import sqrt
from Kspace import Kspace

i = sqrt(-1)
H = Hadamardgate()
I = I()
v0 = basis(2,0)
v1 = basis(2,1)
K=Kspace(9)

def R(G,t):
    return cos(t/2.)*I - i*sin(t/2)*G

def Bounds(A,G,t):
    U = Ugate(A)
    U = Qobj(U)
    B11 = tensor(I,I,I)
    B12 = tensor(I,H,I)
    B23 = tensor(I,v0*v0.dag(),I)+ tensor(I,v1*v1.dag(),U)
    B34 = tensor(I,H,I)
    B45 = tensor(I,v0*v0.dag(),I)+tensor(R(G,t),v1*v1.dag(),I)
    B56 = tensor(I,H,I)
    B67 = tensor(I,v0*v0.dag(),I)+tensor(I,v1*v1.dag(),U.dag())
    B78 = tensor(I,H,I)
    B89 = tensor(v1*v1.dag(),I,I)
    B99 = tensor(I,I,I)
    return [B11,B12,B23,B34,B45,B56,B67,B78,B89,B99]

def ForwardPropagation(w,A,G,t):
    w = sqrt(w)
    L = []
    for i in range(1,10):
        L.append(w*Bounds(A,G,t)[i])
    return L

def ReversePropagation(l,A,G,t):
    l = sqrt(l)
    L = []
    for i in range(0,9):
        L.append(l*Bounds(A,G,t)[i])
    return L

def SupForward(w,A,G,t):
    L = []
    fp = ForwardPropagation(w,A,G,t)
    for i in range(8):
        L.append(tensor(fp[i],K[i+1]*K[i].dag()))
    L.append(tensor(fp[len(fp)-1],K[8]*K[8].dag()))
    return L

def SupReverse(l,A,G,t):
    L = []
    rp = ReversePropagation(l,A,G,t)
    L.append(tensor(rp[0],K[0]*K[0].dag()))
    for i in range(1,9):
        L.append(tensor(rp[i],K[i-1]*K[i].dag()))
    return L 
           