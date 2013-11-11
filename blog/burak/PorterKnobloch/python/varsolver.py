#
# invpolsolver.py
#
"""
Use odeint to solve differential equations defined by velocityvar in twomode.py
"""

from scipy.integrate import odeint
import twomode
import numpy as np

#Callable function version:
def integrate(x0, p, t, abserror=1.0e-8, relerror=1.0e-6):
	"""
	Takes the initial condition, parameters and the time interval
	returns the result as a series in time.
	"""
	xsol = odeint(twomode.velocityvar, x0, t, args=(p,), atol = abserror, rtol = relerror)
	
	return xsol

if __name__ == "__main__":
	
	#Load parameters:
	p = np.loadtxt('data/parameters.dat')

	#Initial conditions:
	x10=-0.310563473904   
	y10=-0.310563473904
	x20=-0.0419063666849
	y20=-0.38981313744
	
	# ODE solver parameters
	abserr = 1.0e-8
	relerr = 1.0e-6
	stoptime = 11.57
	numpoints = 1157

	# Create the time samples for the output of the ODE solver:
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	
	# Pack up the initial conditions:
	x0=np.zeros(4+4*4)
	x0[0:4] = [x10,y10,x20,y20]
	x0[4:4+4*4] = np.identity(4).reshape(16) # Initial condition for the Jacobian is identity

	# Call the ODE solver
	xsol = odeint(twomode.velocityvar, x0, t, args=(p,), atol = abserr, rtol = relerr)

	#Print the solution
	np.savetxt('data/varsolution.dat', xsol)
	np.savetxt('data/vartime.dat', t)
