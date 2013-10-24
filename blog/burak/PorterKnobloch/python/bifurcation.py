#
# bifurcation.py
#
"""
Generates bifurcation diagram by varying a1 and computing the Poincare
section iteratively
"""
import numpy as np
import onslicesolver
import twomode
import psectgs

#Parameter values:
mu1 = 1 
b1 = -1 
c1 = 1 
mu2 = -1 
a2 = 0 
b2 = 0 
c2 = -1 
e2 = 0

xhatphi0 = [-0.713404949283, -0.713404949283, 0.084155488181, -1.5962565603, 0]
#xhatphi0 = [0.01, 0.01, 0, 0, 0]
rangea1 = np.linspace(0.45, 0.49, num=1000)

#dummy assignment for bif:
bif = np.array([[0,0,0,0], [0,0,0,0]], float)

for a1 in rangea1:

	p = [mu1, a1, b1, c1, mu2, a2, b2, c2, e2]
	
	stoptime = 1000
	numpoints = 50000
	
	# Create the time samples for the output of the ODE solver:
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	#Integrate
	xhatphi=onslicesolver.integrate(xhatphi0, p, t)
	#Discard the first 1000 time units and phi:
	xhat = xhatphi[numpoints/2:numpoints,0:4]
	#Transform to the Gram-Schmidt basis:
	xgs = twomode.ssp2gramschmidt(xhat.transpose())
	t = np.array([t[numpoints/2:numpoints]],float)
	xgst = np.concatenate((t,xgs))
	xgst = xgst.transpose()
	
	ps = psectgs.computeps(xgst, p)
	adummy = np.ones(len(ps))*a1
	adummy = np.array([adummy])
	bif = np.append(bif, np.concatenate((adummy.T,ps[:,1:4]), axis=1),0)

bif = bif[2:len(bif),:]

for biff in bif:
	print biff[0], biff[1], biff[2], biff[3]
