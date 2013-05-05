import numpy as np
from qutip import *
from UnitaryGates import *
from Operators import *
from Kspace import *
from cmath import sqrt

i = sqrt(-1)
K = Kspace(9)
Y = Ygate()
P9 = Projection(9,8)
w = 1
l = 1 - w

##############################################################
v0 = basis(2,0)
v1 = basis(2,1)

"""
    The A operator
"""
A = matrix([[1,2*i],[-2*i,1]])
B = Qobj(A)
b0 = B.eigenenergies()[0]
b1 = B.eigenenergies()[1]

###############################################################
#  the b vector
b = 2*B.eigenstates()[1][0]  + 7*B.eigenstates()[1][1]
b = b/b.norm()
###########################################################
if b0 != 0 and b1 != 0:
    if np.abs(b0)>np.abs(b1):
        theta = -2*np.arccos(b1/b0)
    if np.abs(b1)>=np.abs(b0):
        theta = -2*np.arccos(b0/b1)
######################################################################

# density matrix     ###########################
p0 = tensor(v1*v1.dag(),v0*v0.dag(),b*b.dag(),K[0]*K[0].dag())
#p0 = tensor(p0,K[0]*K[0].dag())
P = [p0]
############ M operators ##################################################
Kf = SupForward(w,A,Y,theta)
Kb = SupReverse(l,A,Y,theta)
exit_status = False
n = 0
######################################################################
while not exit_status:
    g = 0
    for j in range(9):
        g += Kf[j]*P[n]*Kf[j].dag() + Kb[j]*P[n]*Kb[j].dag()
    P.append(g)
    n = n + 1
    
    if n>2 and P[len(P)-2]==P[len(P)-1]:
        exit_status = True
    else:
            exit_status = False
 ##########################################################
solution = ptrace(P[len(P)-1]*P9,2)
print solution,n


