#
# movingframes.py
#
"""
Compute the moving frames projection for the full statespace solution 
stored in ../data/sspsolution.dat
"""
import numpy as np
from scipy.optimize import fsolve
import twomode

#Callable function version:

xhatp = np.array([1,0,0,0],float)
T = twomode.generator()
tp = np.dot(T, xhatp)
		
nofphi = 20
deltaphi = 2*np.pi/(nofphi-1)

def SliceCondition(phi, x):
	U = np.dot(tp, np.dot(twomode.LieElement(-phi), x));
	return U

def x2xhatphi(x):
	
	if np.size(x) == 4:
		
		x = x.reshape(1,-1)
	
	xhatphi = np.zeros((np.size(x,0),5), float)
	
	for i in range(len(x)):
		for phi in np.linspace(0, 2*np.pi, num=nofphi):
			if SliceCondition(phi, x[i,:]) > 0 and SliceCondition(phi+deltaphi, x[i,:])<0:
				xhatphi[i,4] = fsolve(SliceCondition, phi, args=(x[i,:],))
				xhatphi[i,0:4] = np.dot(twomode.LieElement(-xhatphi[i,4]),x[i,:]);
				break
	
	return xhatphi

if __name__ == "__main__":
	
	t, x1, x2, y1, y2  = np.loadtxt('data/sspsolution.dat', unpack=True)
	x = np.array([x1,x2,y1,y2], float)
	x=x.transpose()
	
	xhatphi = x2xhatphi(x)
							
	#Print the solution
	for t1, x1 in zip(t,xhatphi):
		print t1, x1[0], x1[1], x1[2], x1[3], x1[4]
	
