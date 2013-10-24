#
#psectgs.py
#
import numpy as np
from scipy.integrate import odeint
import twomode
import onslicesolver

#Load parameters:
pars = np.loadtxt('data/parameters.dat')

xhatpsect = np.array([0,0,0], float)
nhat = np.array([0,1,0], float)
direction = 1 # set 1 or -1 to choose between directions of piercing the Poincare section
tolerance = 1e-9

#Poincare hyperplane equation:
def U(x):
	cond=np.dot(x-xhatpsect, nhat) * direction
	return cond

# ODE solver parameters:
abserr = 1.0e-12
relerr = 1.0e-9

def computeps(xhatGS, p=pars):
	
	deltat = xhatGS[1,0]-xhatGS[0,0] # Read deltat from the data
	
	#dummy assignment for ps:
	ps = np.array([[0, 0, 0, 0], [0,0,0,0]], float)
	
	#Look for intersections:
	for i in range(len(xhatGS)-1):
		if U(xhatGS[i, 1:4])<0 and U(xhatGS[i+1, 1:4])>0:
			xdummyGS = xhatGS[i, 1:4]
			deltatadaptive = deltat/2
			pstime = xhatGS[i,0]
			j=1
			condition = U(xdummyGS)
			xhatphidummy=np.concatenate((twomode.gramschmidt2ssp(xdummyGS), np.array([0],float)))
			while np.abs(condition)>tolerance:
				j=j+1
				#Adaptive integration:
				t = [float(0), deltatadaptive]
				xhatsol = onslicesolver.integrate(xhatphidummy, p, t, abserr, relerr)
				xhatphidummy=xhatsol[1, :]
				pstime = pstime + deltatadaptive
				xdummyGS = twomode.ssp2gramschmidt(xhatphidummy[0:4])
				condition = U(xdummyGS)
				
				if condition > 0:
					xhatphidummy=xhatsol[0,:]
					pstime = pstime - deltatadaptive
					deltatadaptive = deltatadaptive/2
				
				if j > 100:
					print 'not converging, exiting...'
					break
			pstime = np.array([pstime], float)
			psdummy = np.concatenate((pstime, xdummyGS))
		
			ps = np.append(ps, [psdummy],0)
	
	ps = ps[2:len(ps),:]

	return ps

#If the module is called as a script, compute the Poincare section from the
#data and print:

if __name__ == "__main__":
	
	xhatGS = np.loadtxt('data/gramschmidt.dat', unpack=True) #Gram-Schmidt data
	xhatGS=xhatGS.transpose() # xhatGS = [t, x1hatGS, y1hatGS, y2hatGS]
	
	ps=computeps(xhatGS)
	
	#Print the poincare section
	for pss in ps:
		print pss[0], pss[1], pss[2], pss[3]
