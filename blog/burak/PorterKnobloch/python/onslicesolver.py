#
# invpolsolver.py
#
"""
Use odeint to solve differential equations defined by vhatvphi on the symmetry
reduced manifold. Can be called as a script or imported as a module.
"""

import numpy as np 
from scipy.integrate import odeint
import twomode

T = twomode.generator()
xhatp = np.array([1,0,0,0],float)
tp = np.dot(T, xhatp)

def vphi(x,t,p):
	"""
    Velocity function for the group parameter
    """
	vel = np.dot(twomode.vfullssp(x,t,p),tp)/np.dot(np.dot(T,x),tp)
	return vel

def vhatvphi(xphi,t,p):
	"""
    Velocity function within the slice
    """
	vel=np.zeros(5)
	x = xphi[0:4]
	phi = xphi[4]
	vel[0:4] = twomode.vfullssp(x,t,p) - vphi(x,t,p)*np.dot(T,x)
	vel[4] = vphi(x,t,p)
	return vel
	
abserr = 1.0e-8
relerr = 1.0e-6

#Callable function version:
def integrate(xphi0, p, t, abserror=1.0e-8, relerror=1.0e-6):
	"""
	Takes the initial condition, parameters and the time interval
	returns the result as a series in time.
	"""
	xsol = odeint(vhatvphi, xphi0, t, args=(p,), atol = abserror, rtol = relerror)
	
	return xsol


#Integrate only if the module is called as a script:

if __name__ == "__main__":
	
	#Load parameters:
	p = np.loadtxt('data/parameters.dat')

	#Initial conditions:
	
	#x10=-0.311102644920502
	#y10=-0.311102644920502
	#x20=-0.07020440068369531
	#y20=-0.3862670595294275
	#phi0= 0
	#Releq:
	x10=0.43997
	y10=0
	x20=-0.38627
	y20=0.07020
	phi0= 0
	
	#RPO:
	#x10=0.4525858904342642
	#y10=0.0000000000000000
	#x20=0.0510590531307565
	#y20=0.0318533099488150

	#RPO:
	#x10=0.4525858904295301
	#y10=0.0000000000000000
	#x20=0.0510590529667102
	#y20=0.0318533099915381
	
	#RPO2:
	#x10=2.22664041e-01  
	#y10=-9.09511319e-08   
	#x20=8.72897436e-02   
	#y20=3.09876107e-02
	
	#RPO2:
	#x10= 0.45169971          
	#y10=0
	#x20= 0.02035061
	#y20=0.03874689
   

	#Xi
	x10=0.45258589   
	y10=0.0 
	x20=0.05105905
	y20=0.03185331
		
	# ODE solver parameters
	abserr = 1.0e-14
	relerr = 1.0e-13
	stoptime = 1000
	stoptime =  3.63
	numpoints = 100000+1
	numpoints = 363
	
	# Create the time samples for the output of the ODE solver:
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	
	# Pack up the initial conditions:
	xphi0 = [x10,y10,x20,y20, phi0]
	
	# Call the ODE solver
	xsol = odeint(vhatvphi, xphi0, t, args=(p,), atol = abserr, rtol = relerr)
	
	#Print the solution
	for t1, x1 in zip(t,xsol):
		print t1, x1[0], x1[1], x1[2], x1[3], x1[4]
