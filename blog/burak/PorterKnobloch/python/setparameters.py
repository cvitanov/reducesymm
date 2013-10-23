#
# setparameters.py
#
import numpy as np

#Parameter values:
mu1 = 1 
a1 = 0.47 
b1 = -1 
c1 = 1 
mu2 = -1 
a2 = 0 
b2 = 0 
c2 = -1 
e2 = 0

p = [mu1, a1, b1, c1, mu2, a2, b2, c2, e2]
np.savetxt('data/parameters.dat', p)
