#
#psectgs.py
#
import numpy as np
from scipy.integrate import odeint
import twomode

#Load parameters:
p = np.loadtxt('data/parameters.dat')

xhatpsect = np.array([0,0,0], float)
nhat = np.array([0,1,0], float)
direction = 1 # set 1 or -1 to choose between directions of piercing the Poincare section
tolerance = 1e-9

xhatGS = np.loadtxt('data/gramschmidt.dat', unpack=True) #Gram-Schmidt data
xhatGS=xhatGS.transpose() # xhatGS = [t, x1hatGS, y1hatGS, y2hatGS]
deltat = xhatGS[1,0]-xhatGS[0,0] # Read deltat from the data
xhat = np.loadtxt('data/solutiononslice.dat', unpack=True) #Reduced ssp data
xhat = xhat.transpose() # xhat = [t, x1hat, x2hat, y1hat, y2hat, phi]

#On slice evolution equations:

T = twomode.generator()
xhatp = np.array([1,1,0,0],float)
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
	
#Poincare hyperplane equation:

def U(x):
	cond=np.dot(x-xhatpsect, nhat) * direction
	return cond

# ODE solver parameters:
abserr = 1.0e-12
relerr = 1.0e-9
#stoptime = 1000
#numpoints = 100000

#dummy assignment for ps:
ps = np.array([[0, 0, 0, 0], [0,0,0,0]], float)

#Look for intersections:
for i in range(len(xhatGS)-1):
	if U(xhatGS[i, 1:4])<0 and U(xhatGS[i+1, 1:4])>0:
		xdummyGS = xhatGS[i, 1:4]
		deltatadaptive = deltat/2
		pstime = i*deltat;
		j=1	;
		condition = U(xdummyGS)
		xhatphidummy=xhat[i, 1:6]
		while np.abs(condition)>tolerance:
			j=j+1
			#Adaptive integration:
			t = [float(0), deltatadaptive]
			xhatsol = odeint(vhatvphi, xhatphidummy, t, args=(p,), atol = abserr, rtol = relerr)
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

#Print the poincare section
for pss in ps:
	print pss[0], pss[1], pss[2], pss[3]
