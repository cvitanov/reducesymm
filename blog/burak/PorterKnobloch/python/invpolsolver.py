#
# invpolsolver.py
#
"""
Use odeint to solve differential equations defined by vinvpol in twomode.py
"""

from scipy.integrate import odeint
import twomode
import numpy as np

#Load parameters:
p = np.loadtxt('data/parameters.dat')

#Initial conditions:
u0 = 1.5652
v0 = 2.5564
w0 = -4.3934
q0 = 1.7731

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 1000
numpoints = 100000

# Create the time samples for the output of the ODE solver:
t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

# Pack up the parameters and initial conditions:
x0 = [u0,v0,w0,q0]

# Call the ODE solver
xsol = odeint(twomode.vinvpol, x0, t, args=(p,), atol = abserr, rtol = relerr)

#Print the solution
for t1, x1 in zip(t,xsol):
	print t1, x1[0], x1[1], x1[2], x1[3]
