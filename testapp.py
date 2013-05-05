import numpy as np
from qutip import *
from UnitaryGates import *
from Operators import *
from Kspace import *
import cmath
##############################################################
i = cmath.sqrt(-1)
K = Kspace(9)
Y = Ygate()
P9 = Projection(9,8)
w = 0.6
l = 1-w

v0 = basis(2,0)
v1 = basis(2,1)
###################################################################
"""
    The A operator
"""
A = matrix([[0.5,0.32*i],[-0.2*i,0.11]])
B = Qobj(A)
b0 = B.eigenenergies()[0]
b1 = B.eigenenergies()[1]

###############################################################
#  the b vector expressed in terms of the eigenvectors of A
u1 = B.eigenstates()[1][0]
u2 = B.eigenstates()[1][1]
b = 0.2*u1  + 0.7*i*u2
b = b/b.norm()
###########################################################
if b0 != 0 and b1 != 0:
    if np.abs(b0)>np.abs(b1):
        theta = -2*np.arccos(b1/b0)
    if np.abs(b1)>=np.abs(b0):
        theta = -2*np.arccos(b0/b1)
else:
    theta = np.pi
######################################################################
"""
    Inital density matrix
"""
p0 = tensor(v1*v1.dag(),v0*v0.dag(),b*b.dag(),K[0]*K[0].dag())
P = [p0]
#####################################################################
"""
    The M_j^i operators defined from the Operators Module
"""
Kf = SupForward(w,A,Y,theta)
Kb = SupReverse(l,A,Y,theta)
exit_status = False
n = 0
######################################################################
while not exit_status:
    g = tensor(Z(2),Z(2),Z(2),Z(9))
    for j in range(9):
        g = g + Kf[j]*P[n]*Kf[j].dag() + Kb[j]*P[n]*Kb[j].dag()
    P.append(g)
    n = n + 1
    
    if n>11 and P[len(P)-2]*P9==P[len(P)-1]*P9:
        exit_status = True
    else:
            exit_status = False
 ##########################################################
steadystate = P[len(P)-1]
print steadystate,n



