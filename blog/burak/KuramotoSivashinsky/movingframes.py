#
# movingframes.py
#
"""
Compute the moving frames projection for the full statespace solution 
stored in ../data/sspsolution.dat
"""
import numpy as np
from scipy.optimize import fsolve
from scipy.linalg import expm
import KS

N=15;
xhatp = np.zeros(2*N,float)
xhatp[0]=1
xhatp[1]=0
T = KS.generator(N)
tp = np.dot(T, xhatp)

xd = np.loadtxt('data/sspsolution.dat')
t = np.zeros(np.size(xd,0))
t = xd[:,0]
x = np.zeros((np.size(xd,0), np.size(xd,1)-1), float)
x = xd[:,1:np.size(xd,1)]

#x=x.transpose()

def LieElement(phi):
	g = expm(phi*T)
	return g
	

def SliceCondition(phi, x):
	U = np.dot(tp, np.dot(LieElement(-phi), x))*(-1);
	return U

nofphi = 20
deltaphi = 2*np.pi/(nofphi-1)

xhatphi = np.zeros((len(x),np.size(x,1)+1), float)
for i in range(len(x)):
	for phi in np.linspace(0, 2*np.pi, num=nofphi):
		if SliceCondition(phi, x[i,:]) < 0 and SliceCondition(phi+deltaphi, x[i,:])>0:
			xhatphi[i,np.size(xhatphi, 1)-1] = fsolve(SliceCondition, phi, args=(x[i,:],))
			xhatphi[i,0:np.size(xhatphi, 1)-1] = np.dot(LieElement(-xhatphi[i,np.size(xhatphi, 1)-1]),x[i,:]);
			break
			
			
#Print the solution
np.savetxt('data/movingframes.dat', xhatphi)
np.savetxt('data/time.dat', t)
#for t1, x1 in zip(t,xhatphi):
#	print t1, x1[0:np.size(x1,0)-1] #, x1[1], x1[2], x1[3], x1[4]
