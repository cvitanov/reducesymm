#
# bifurcation.py
#
"""
Generates bifurcation diagram by varying a1 and computing the Poincare
section iteratively
"""
import numpy as np
import invpolsolver
import twomode

#Parameter values:
mu1 = 1 
b1 = -1 
c1 = 1 
e1 = 0
mu2 = -1 
a2 = 0 
b2 = 0 
c2 = -1 
e2 = 0

u0 = [2, 1, 1, 3.872983]

rangea1 = np.linspace(0.45, 0.485, num=100)

#dummy assignment for bif:
bif = np.array([[0,0], [0,0]], float)

for a1 in rangea1:

	p = [mu1, a1, b1, c1, e1, mu2, a2, b2, c2, e2]
	
	stoptime = 1000
	numpoints = 50000
	
	# Create the time samples for the output of the ODE solver:
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	#Integrate
	ut=invpolsolver.integrate(u0, p, t)
	#Discard the first 1000 time units:
	#utsteady = ut[numpoints/2:numpoints,0:4]
	qtsteady = np.copy(ut[numpoints/10:numpoints,3])
	qtsteady = np.array([qtsteady])
	
	adummy = np.ones(np.size(qtsteady,1))*a1
	adummy = np.array([adummy])
	bif = np.append(bif, np.concatenate((adummy.T,qtsteady.T), axis=1),0)

bif = bif[2:len(bif),:]

for biff in bif:
	print biff[0], biff[1]
