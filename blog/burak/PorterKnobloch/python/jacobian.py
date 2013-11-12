#
# jacobian.py
#
"""
Loads the solution of the variational equation for one period.
Computes the Jacobian and eigenvalues for the period.
"""
import numpy as np

x = np.loadtxt('data/varsolution.dat')

J = x[np.size(x,0)-1, 4:20].copy().reshape(4,4)
w,v = np.linalg.eig(J)

print w
print v
