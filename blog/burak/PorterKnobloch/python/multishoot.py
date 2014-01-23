#
#multishoot.py
#
"""
Multiple shooting solver to find relative periodic orbits
"""

import numpy as np
from scipy import interpolate, integrate
from scipy.misc import derivative
from scipy.optimize import newton, fsolve, root
import numdifftools as nd

import twomode
import onslicesolver
import sspsolver

#Read Porter - Knobloch system parameters:

pars = np.loadtxt('data/parameters.dat')

#Load the RPO candidates:

group = np.loadtxt('data/group.dat')
itineraries = np.loadtxt('data/itineraries.dat', dtype="str")
position = np.loadtxt('data/position.dat')
tof0 = np.loadtxt('data/tof0.dat')
phi0 = np.loadtxt('data/phi0.dat')
xrpo0 = np.loadtxt('data/xrpo0.dat')

#Let's handle the first one:
#Compute normal vector to the Poincare section:
sectp = np.array([0.4399655797367152, 0, -0.38626706847930564, 0.0702043939917171], float)
sectp3D = np.array([sectp[0], sectp[2], sectp[3]],float)
vaux = np.array([0,0,1], float) #Auxiliary vector to decide the Poincare section direction
unstabledir = np.array([0.02884567,  0.99957642, -0.00386173], float) #Unstable direction
nhat3D = np.cross(vaux, unstabledir)
nhat = np.array([nhat3D[0],0,nhat3D[1],nhat3D[2]], float)

#Initial guess:
x0 = xrpo0[0,:]
T0 = tof0[0]

xt0 = np.append(x0, tof0[0])

xhatphi = np.append(x0, np.array([0], float))

def ftau(xhat, tau):
	"""
	Lagrangian description of the reduced flow.
	"""
	#Add the group parameter to be able to integrate:
	xhatphi = np.append(xhat, np.array([0], float))
	stoptime = tau
	numpoints = 2
	
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	
	xsol = onslicesolver.integrate(xhatphi, pars, t, abserror=1.0e-12, relerror=1.0e-11)
	xhattau = xsol[1,0:4]
	
	return xhattau

converged = False
tol = 1e-9
xi = x0
Ti = T0

#Define maximum number of steps:
imax = 10
#Apply ChaosBook p294 (13.11) 
#with constraint nhat . Dx = 0
i=0
tol = 1e-9

#while not converged:
	
fx = lambda x: ftau(x, Ti)

xTi = ftau(xi, Ti)
error = xi - xTi
print "error: "
print error

error0 = np.append(error, np.array([0]))

if max(np.abs(error)) < tol:
	
	xrpo = xi
	Trpo = Ti
	converged = True
	#break

Jfx = nd.Jacobian(fx)
J = Jfx(xi)

print "Jacobian:"
print J

vx = onslicesolver.vhatvphi(np.append(xTi, np.array([0], float)),0,pars)[0:4]
vx = vx.reshape(-1,1)

print "vx:"
print vx

A = np.concatenate((np.identity(4)-J, vx), axis=1)

A = np.concatenate((A, np.append(nhat.reshape(-1,1), np.array([[0]], float), axis=0).transpose()), axis=0)

np.savetxt('Amatrix', A)

Ainv = np.linalg.inv(A)
print "Ainv:"
print Ainv

print "A.Ainv:"
print np.dot(A, Ainv)

DxT = np.dot(Ainv, -error0)

print "DxT"
print DxT

xi = xi + DxT[0:4]
Ti = Ti + DxT[4]

xTi = ftau(xi, Ti)
error = xi - xTi
print "error: "
print error

i = i + 1

	#if i > imax:
		
		#print "did not converged in given maximum number of steps"
		#print "exitting..."
		#xrpo = xi
		#Trpo = Ti
		#break
	
	#print error
	
print "x_rpo:"
print xrpo
print "T_rpo:"
print Trpo
print "Error"
print error
