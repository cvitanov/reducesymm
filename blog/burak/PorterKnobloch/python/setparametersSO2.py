#
# setparameters.py
#
import numpy as np

#Parameter values:
mu1 = -2.8 
a1 = -1
b1 = 0
c1 = -7.75
e1 = 0
mu2 = 1 
a2 = -2.66 
b2 = 0 
c2 = 1
e2 = 1

p = [mu1, a1, b1, c1, e1, mu2, a2, b2, c2, e2]
np.savetxt('data/parameters.dat', p)
