#
# setparameters.py
#
import numpy as np

#Parameter values:
mu1 = 1 
a1 = 0.452121212121
b1 = -1 
c1 = 1 
e1 = 0
mu2 = -1 
a2 = 0 
b2 = 0 
c2 = -1 
e2 = 0

#p = [mu1, a1, b1, c1, mu2, a2, b2, c2, e2]
p = [mu1, a1, b1, c1, e1, mu2, a2, b2, c2, e2]
np.savetxt('data/parameters.dat', p)
