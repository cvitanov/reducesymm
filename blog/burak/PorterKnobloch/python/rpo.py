#
#INCOMPLETE. SHOULD NOT BE USED AND LIKELY TO PRODUCE WRONG RESULTS. 
#HOWEVER, CAN INCLUDE APPROACHES THAT MAY BE USEFUL FOR FUTURE APPLICATIONS
#
#psectandretmap.py
#
"""
	Calculates the poincare section intersections using the previously computed
	solution on the slice.
	Plots the poincare section.
	Computes and plots the return map.
	This will produce correct results only if the xhatp=[1,0,0,0] slice 
	is used.
"""

import numpy as np
from scipy.optimize import fsolve, newton_krylov, anderson

import numpy as np

import twomode
import onslicesolver
#import psectslice

#Load parameters:
p = np.loadtxt('data/parameters.dat')

def rpo(x, Tapproximate):
	print "x"
	print x
	abserr = 1.0e-14
	relerr = 1.0e-12
	starttime = Tapproximate*0.9
	stoptime = Tapproximate*1.1
	#print "starttime = ", starttime
	#print "stoptime = ", stoptime
	#numpoints = int(stoptime/0.0001 + 1.0)
	numpoints1 = int(starttime/0.01 + 1.0)
	numpoints = int(1e4)
	#print "numpoints:"
	#print numpoints
	#numpoints = 2
	#xphi0=[x[0], x[1], x[2], x[3], 0]
	xphi0=[x[0], 0, x[1], x[2], 0]
	
	t1=[starttime * float(i) / (numpoints1 - 1) for i in range(numpoints1)]
	xsol1 = onslicesolver.integrate(xphi0, p, t1, abserror=abserr, relerror=relerr)
	xphi1 = xsol1[len(xsol1)-1,:]
	
	# Create the time samples for the output of the ODE solver:
	#t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	#t = np.append([0], np.linspace(starttime, stoptime, numpoints))
	
	t = [(stoptime-starttime) * float(i) / (numpoints - 1) for i in range(numpoints)]
	#xsol = onslicesolver.integrate(xphi0, p, t, abserror=abserr, relerror=relerr)
	xsol = onslicesolver.integrate(xphi1, p, t, abserror=abserr, relerror=relerr)
	xsol3d = np.concatenate((np.array([xsol[:,0]]).T, xsol[:,2:4]), axis=1)
	
	#diff = xsol[int(numpoints/2):len(xsol),0:4] - x
	#diff = xsol3d[int(numpoints/2):len(xsol),0:3] - x
	diff = xsol3d[0:len(xsol3d),0:3] - x
	
	diff = np.linalg.norm(diff, axis=1)
	
	#imin = np.argmin(diff)+int(numpoints/2)
	#imin = np.argmin(diff)+1
	imin = np.argmin(diff)
	
	print "tmin"
	print t[imin]+starttime
	
	
	
	#xf = [xsol[imin, 0], 
		  #xsol[imin, 1], 
		  #xsol[imin, 2], 
		  #xsol[imin, 3]] 
	xf = [xsol3d[imin, 0], 
		  xsol3d[imin, 1], 
		  xsol3d[imin, 2]] 
	
	print "xf"
	print xf
	
	return xf - x
	
def findrpo(x0, Tapproximate):
	rpos = fsolve(rpo, x0, args=(Tapproximate,))
	return rpos
	
def rpot(xt):
	'''
	not for use, under development
	'''
	#print "xt"
	#print xt
	abserr = 1.0e-14
	relerr = 1.0e-12
	#stoptime = xt[4]
	stoptime = xt[3]
	numpoints = int(stoptime/0.01 + 1.0)
	#print "numpoints:"
	#print numpoints
	#numpoints = 2
	#xphi0=[xt[0], xt[1], xt[2], xt[3], 0]
	xphi0=[xt[0], 0, xt[1], xt[2], 0]
	
	# Create the time samples for the output of the ODE solver:
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	
	xsol = onslicesolver.integrate(xphi0, p, t, abserror=abserr, relerror=relerr)
	xsol3d = np.concatenate((np.array([xsol[:,0]]).T, xsol[:,2:4]), axis=1)
	
	#diff = np.concatenate((xsol[10:len(xsol),0:4], np.array([t[10:len(xsol)]]).T), axis=1) - xt[0:5]
	diff = np.concatenate((xsol3d[10:len(xsol),0:3], np.array([t[10:len(xsol)]]).T), axis=1) - xt[0:5]
	#diff = np.linalg.norm(diff, axis=1)

	#diff = xsol[10:len(xsol),0:4] - xt[0:4]
		
	#diff[:,4] = diff[:,4]/xt[4]
	diff = np.linalg.norm(diff, axis=1)
	
	imin = np.argmin(diff)+10
	
	#print "imin"
	#print imin
	
	#xf = [xsol[imin, 0], 
		  #xsol[imin, 1], 
		  #xsol[imin, 2], 
		  #xsol[imin, 3], 
		  #t[imin]] 
	xf = [xsol3d[imin, 0], 
		  xsol3d[imin, 1], 
		  xsol3d[imin, 2], 
		  t[imin]] 
	
	#print "xf"
	#print xf
	
	return xf - xt

def findrpoT(xt0):
	'not ready for use'
	rpos = fsolve(rpot, xt0)
	return rpos

def findrpoNK(x0, Tapproximate):
	
	def rpot(x):
		return rpo(x, Tapproximate)
		
	
	rpos = newton_krylov(rpot, x0)
	return rpos

def findrpoA(x0, Tapproximate):
	
	def rpot(x):
		return rpo(x, Tapproximate)
		
	
	rpos = anderson(rpot, x0)
	return rpos
	
if __name__ == "__main__":

	#candidate:
	x10=0.4525858904295301
	y10=0.0000000000000000
	x20=0.0510590529667102
	y20=0.0318533099915381
	phi0=0
	t0 = 3.63
	
	xt0 = np.array([x10, y10, x20, y20, t0], float)
	
	xt0 = np.array([4.551491852618429479e-01, 0.000000000000000000e+00, 1.398837906790846852e-01, 4.074262645591256171e-03, 1.096855715200001669e+01], float)
	xt0 = np.array([4.525769625920411654e-01, 0.000000000000000000e+00, 5.074968050853578827e-02, 3.193047332800166838e-02, 3.650669934999996258e+00], float)
	xt03d = np.array([4.525769625920411654e-01, 5.074968050853578827e-02, 3.193047332800166838e-02, 3.650669934999996258e+00], float)
	xt03d = np.array([4.537403786922851179e-01, 9.106503114702829693e-02, 2.080252458763354878e-02, 2.204607994399998461e+01], float)
	xt03d = np.array([4.533343012015041973e-01, 7.699340499040213670e-02, 2.494053433095596434e-02, 7.347242386999937480e+00], float)
	xt03d = np.array([4.551491852618429479e-01, 1.398837906790846852e-01, 4.074262645591256171e-03, 1.096855715200001669e+01], float)
	#xt03d = np.array([4.535000042007177878e-01, 8.273543874233757478e-02, 2.328624944931542706e-02, 2.588854215399999248e+01], float)
	
	#rpos = findrpo(xt0[0:4], xt0[4]) # 4D
	rpos = findrpo(xt03d[0:3], xt03d[3]) # 3D Newton
	#rpos = findrpoNK(xt03d[0:3], xt03d[3]) # 3D Newton - Krylov
	#rpos = findrpoA(xt03d[0:3], xt03d[3]) # 3D Anderson
	
	print "rpo at"
	print rpos
	print "fval:"
	print rpo(rpos, xt03d[3])
