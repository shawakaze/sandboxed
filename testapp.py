import numpy as np
from qutip import *
from UnitaryGates import *
from Operators import *
from Kspace import *

K = Kspace(9)
Y = Ygate()
w = 0.8
l = 1 - w
pi = np.pi

##############################################################
v0 = basis(2,0)
v1 = basis(2,1)
b = v1
# density matrix     ############################
p0 = tensor(v1*v1.dag(),v0*v0.dag(),b*b.dag())
P = [p0]
############## Operators ########################################################
Bf = ForwardPropagation(w,Y,pi/2.)
Br = ReversePropagation(l,Y,pi/2.)
##################################################################################

for j in range(9):
    g = Bf[j]*P[j]*Bf[j].dag() + Br[j]*P[j]*Br[j].dag()
    P.append(g)


solution = P[len(P)-1]

