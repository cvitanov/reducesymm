#
#psectslice.py
#
import numpy as np
from scipy.integrate import odeint
import twomode

#Load parameters:
pars = np.loadtxt('data/parameters.dat')

# ODE solver parameters:
abserr = 1.0e-14
relerr = 1.0e-13
#Zero olerance:
tolerance = 1e-14

def computeps(xsol, sectp, nhat, direction,  p=pars):
	
	#Poincare hyperplane equation:
	def U(x):
		cond=np.dot(x-sectp, nhat) * direction
		return cond
	
	#Apply hyperplane condition to each element:
	uxarray = np.dot(xsol[:,1:5] - sectp, nhat) * direction
	
	deltat = xsol[1,0]-xsol[0,0] # Read deltat from the data
	
	#dummy assignment for ps:
	ps = np.array([[0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], float)
	
	#Look for intersections:
	#for i in range(len(uxarray)-1):
		#if uxarray[i]<0 and uxarray[i+1]>0:
	for i in range(len(xsol)-1):
		if U(xsol[i, 1:5])<0 and U(xsol[i+1, 1:5])>0:	
			
			xhatphidummy = xsol[i, 1:6]
			deltatadaptive = deltat/2
			pstime = xsol[i,0]
			j=1
			condition = U(xhatphidummy[0:4])
			while np.abs(condition)>tolerance:
				j=j+1
				#Adaptive integration:
				t = [float(0), deltatadaptive]
				xhatsol = twomode.intslice(xhatphidummy, t, p, abserr, relerr)
				xhatphidummy=xhatsol[1, :]
				pstime = pstime + deltatadaptive
				
				condition = U(xhatphidummy[0:4])
				
				if condition > 0:
					xhatphidummy=xhatsol[0,:]
					pstime = pstime - deltatadaptive
					deltatadaptive = deltatadaptive/2
				
				if j > 100:
					print 'not converging, exiting...'
					break
			pstime = np.array([pstime], float)
			psdummy = np.concatenate((pstime, xhatphidummy[0:4]))
		
			ps = np.append(ps, [psdummy],0)
	
	ps = ps[2:len(ps),:]

	return ps

#If the module is called as a script, compute the Poincare section from the
#data and print:

if __name__ == "__main__":
		
	sectp = np.array([0.4399655797367152, 0, -0.38626706847930564, 0.0702043939917171], float)
	vdummy = np.array([0,0,1], float)
	unstabledir = np.array([0.02884567,  0.99957642, -0.00386173], float)
	nhatdummy = np.cross(vdummy, unstabledir)
	
	nhat = np.array([nhatdummy[0],0,nhatdummy[1],nhatdummy[2]], float)
	direction = 1 # set 1 or -1 to choose between directions of piercing the Poincare section
	
	xhat = np.loadtxt('data/solutiononslice.dat') #Solution on slice
	
	ps=computeps(xhat, sectp, nhat, direction,  p=pars)
	
	#Print the poincare section
	for pss in ps:
		print pss[0], pss[1], pss[2], pss[3], pss[4]
