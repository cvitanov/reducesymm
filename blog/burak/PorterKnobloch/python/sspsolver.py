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
#x10 = 0.181612
#x20 = -0.512500
#y10 = 1.939339
#y20 = -0.010343
x10 = 0.01
x20 = 0.01
y10 = 0.01
y20 = 0.01

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 1000
numpoints = 100000

# Create the time samples for the output of the ODE solver:
t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

# Pack up the initial conditions:
x0 = [x10,x20,y10,y20]

# Call the ODE solver
xsol = odeint(twomode.vfullssp, x0, t, args=(p,), atol = abserr, rtol = relerr)

#Print the solution
for t1, x1 in zip(t,xsol):
	print t1, x1[0], x1[1], x1[2], x1[3]
