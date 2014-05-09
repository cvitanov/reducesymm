#
#multishoot.py
#
"""
Multiple shooting solver to find relative periodic orbits
"""
import sys

import numpy as np
from scipy import interpolate, integrate
from scipy.misc import derivative
from scipy.optimize import newton, fsolve
import numdifftools as nd

import twomode
import sspsolver
import varsolver

from pylab import figure, grid, hold, legend, savefig, plot, xlabel, ylabel, show

#Read Porter - Knobloch system parameters:

pars = np.loadtxt('data/parameters.dat')

#Load the RPO candidates:

group = np.loadtxt('data/group.dat')
itineraries = np.loadtxt('data/itineraries.dat', dtype="str")
position = np.loadtxt('data/position.dat')
tof0 = np.loadtxt('data/tof0.dat')
phirpo0 = np.loadtxt('data/phi0.dat')
xrpo0 = np.loadtxt('data/xrpo0.dat')

#Let's handle the first one:
#Compute normal vector to the Poincare section:
sectp = np.array([0.4399655797367152, 0, -0.38626706847930564, 0.0702043939917171], float)
sectp3D = np.array([sectp[0], sectp[2], sectp[3]],float)
vaux = np.array([0,0,1], float) #Auxiliary vector to decide the Poincare section direction
unstabledir = np.array([0.02884567,  0.99957642, -0.00386173], float) #Unstable direction
nhat3D = np.cross(vaux, unstabledir)
nhat = np.array([nhat3D[0],0,nhat3D[1],nhat3D[2]], float)

#The slice:
T = twomode.generator()
xhatp = np.array([1,0,0,0],float)
tp = np.dot(T, xhatp)

#Initial guess:
x0 = xrpo0[0,:]
T0 = tof0[0]
phi0 = -phirpo0[0]

xt0 = np.append(x0, tof0[0])
xtphi0 = np.append(xt0, -phirpo0[0])

def ftau(x, tau):
	"""
	Lagrangian description of the reduced flow.
	"""

	stoptime = tau
	numpoints = 2
	
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	
	xsol = sspsolver.integrate(x, pars, t, abserror=1.0e-10, relerror=1.0e-9)
	fxtau = xsol[1,:]
	
	return fxtau

def Jacobian(x, Ti):
	"""
	Jacobian J^Ti(x)
	"""
	#n = np.size(x,0)
	#J = np.zeros((n,n))
	#for i in range(n):
		#for j in range(n):
			#def fxji(xj):
				
				#xx = x
				#xx[j] = xj
				#fxji = ftau(xx, Ti)[i]
				#return fxji
			
			#J[i,j] = derivative(fxji, x[j])
	
	xvar = np.append(x, np.identity(4).reshape(16))

	stoptime = Ti
	numpoints = 2
	
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	
	xvarsol = varsolver.integrate(xvar, pars, t, abserror=1.0e-10, relerror=1.0e-9)

	J = xvarsol[1, 4:20].reshape(4,4)

	return J		

tol = 1e-9
earray = np.array([], float)
xrpo = np.zeros(np.shape(xrpo0))
tofrpo = np.zeros(np.shape(tof0))
phirpo = np.zeros(np.shape(phirpo0))

Adaptive = True
iAdaptiveMax = 10			
factor = 2
			
for i in range(1, int(np.max(group)+1)):

	AdaptiveFail = False
	print "Group no: ", i
	
	#Get corresponding indices for the group i:
	gindices = np.argwhere(group == i)
	#Reshape them in a 1D array:
	gindices = gindices.reshape(np.size(gindices))
	#Get coordinates of group i
	xi = xrpo0[gindices,:]
	#Get time of flights for group i:
	Ti = tof0[gindices]
	#Get group parameter for group i:
	phii = -phirpo0[gindices]
	#Get the positions in group i
	posi = np.array([position[gindices]-1], int)
	posi = posi.reshape(np.size(posi))
	print "posi = ", posi 
	#Sort with respect to positions:
	xi = xi[posi, :]
	print "xi = ", xi 
	Ti = np.array([Ti[posi]])
	Ti = Ti.reshape(np.size(Ti))
	print "Ti = ", Ti 
	phii = np.array([phii[posi]])
	phii = phii.reshape(np.size(phii))%(2*np.pi)
	print "phii = ", phii 
	
	#raw_input("Press Enter to continue...")
	
	#How many points:
	npts = np.size(gindices)
	#Num of dim:
	n = 4
	#Initiate A, xT, xTT, error, error00 matrices:
	A = np.zeros(((n+2)*npts, (n+2)*npts))
	xT = np.zeros((npts, n))

	error = np.zeros((npts, n))
	error00 = np.zeros(npts*(n+2))
	
	#Dummy variables used in the adaptive block:
	xTT = np.zeros(np.shape(xT))
	xii = np.zeros(np.shape(xi))
	Tii = np.zeros(np.shape(Ti))
	phiii = np.zeros(np.shape(Ti))
	errorr = np.zeros(np.shape(error))
	
	#Define maximum number of iterations:
	itermax = 200
	#Apply ChaosBook p294 (13.11) 
	#with constraint nhat . Dx = 0
	iteration=0
	converged = False
	earray = np.array([], float)
	
	while not(converged):
		
		print "xi"
		print xi
		#print "phii"
		#print phii
		#print "Ti"
		#print Ti
		
		A = np.zeros(((n+2)*npts, (n+2)*npts))
		for j in range(npts):
			
			xTj = ftau(xi[j], Ti[j])
			xT[j, :] = xTj
			
			#Construct Aj matrix for each cycle element:
			#Jacobian:
			J = Jacobian(xi[j], Ti[j])
			
			#minusgvx term:
			vx = np.array(twomode.vfullssp(xTj,0,pars), float)
			minusgvx = -np.dot(twomode.LieElement(phii[j]), vx)
			minusgvx = minusgvx.reshape(-1,1)
			
			#minusTgfTx term:
			minusTgfTx = np.dot(twomode.LieElement(phii[j]), xTj)
			minusTgfTx = -np.dot(T, minusTgfTx)
			minusTgfTx = minusTgfTx.reshape(-1,1)
			
			
			Aj = np.concatenate((-np.dot(twomode.LieElement(phii[j]), J), minusgvx), axis=1)
			Aj = np.concatenate((Aj, minusTgfTx), axis=1)
			
			#nhat0 = np.append(nhat.reshape(-1,1), np.array([[0], [0]], float))
			
			#print nhat0
			
			Aj = np.concatenate((Aj, np.append(vx.reshape(-1,1), np.array([[0], [0]], float), axis=0).transpose()), axis=0)
			#Aj = np.concatenate((Aj, np.append(nhat.reshape(-1,1), np.array([[0], [0]], float), axis=0).transpose()), axis=0)
			Aj = np.concatenate((Aj, np.append(tp.reshape(-1,1), np.array([[0], [0]], float), axis=0).transpose()), axis=0)		
			
			#print Aj

			A[j*(n+2):(j+1)*(n+2), j*(n+2):(j+1)*(n+2)] = Aj
		
			A[j*(n+2):j*(n+2)+n, ((j+1)%npts)*(n+2) : ((j+1)%npts)*(n+2) + n] = A[j*(n+2):j*(n+2)+n, ((j+1)%npts)*(n+2) : ((j+1)%npts)*(n+2) + n] + np.identity(n)
		
		for k in range(npts):
			for l in range(npts):
				print "A("+str(k)+", "+str(l)+")"
				print A[k*(n+2):(k+1)*(n+2), l*(n+2):(l+1)*(n+2)]
					
		for j in range(npts):
			
			error[j, :] = xi[j] - np.dot(twomode.LieElement(phii[(j-1)%npts]), xT[(j-1)%npts,:])
			error00[((j-1)%npts)*(n+2):((j-1)%npts + 1)*(n+2)] = np.append(error[j,:], np.array([0, 0]))
		
		print "error:"
		print error
		
		earray = np.append(earray, np.max(np.abs(error)))
					
		Ainv = np.linalg.inv(A)
		
		print "Ainv.A"
		print np.dot(Ainv, A)
		#idntty = np.dot(Ainv, A)
	
		DxTphi = np.dot(Ainv, -error00)
		
		alpha = 1
		iadaptive = 0
		
		if Adaptive:
			
			Converging = False
			
			while not(Converging):
				
				DxTphi = DxTphi*alpha
				iadaptive += 1
				print "iadaptive= ", iadaptive
				print "alpha= ", alpha
				
				for j in range(npts):
			
					xii[j,:] = xi[j,:] + DxTphi[j*(n+2):j*(n+2)+n]
					Tii[j] = Ti[j] + DxTphi[j*(n+2)+n]
					phiii[j] = phii[j] + DxTphi[j*(n+2)+n+1]			
				
					xTT[j, :] = ftau(xii[j], Tii[j])
	
				for j in range(npts):
					
					errorr[j, :] = xii[j] - np.dot(twomode.LieElement(phiii[(j-1)%npts]), xTT[(j-1)%npts,:])
				
				emax = np.max(np.abs(errorr))
				emaxprevious = earray[len(earray)-1]
					
				if  emax < emaxprevious:
					
					Converging = True
				
				elif iadaptive > iAdaptiveMax:
					
					print "Adaptive step failed, looks kinda hopeless."
					AdaptiveFail = True
					break
				
				else:					
					
					alpha = float(alpha) / float(factor)
					
				#raw_input("Press Enter to continue...")	

		for j in range(npts):
			#print "j ="	
			#print j	
			#print "((j-1)%npts)*(n+2) ="
			#print ((j-1)%npts)*(n+2)
			xi[j,:] = xi[j,:] + DxTphi[j*(n+2):j*(n+2)+n]
			#print "((j-1)%npts)*(n+2)+n ="
			#print ((j-1)%npts)*(n+2)+n
			Ti[j] = Ti[j] + DxTphi[j*(n+2)+n]
			#print "((j-1)%npts)*(n+2)+n+1 ="
			#print ((j-1)%npts)*(n+2)+n+1
			phii[j] = phii[j] + DxTphi[j*(n+2)+n+1]


		iteration += 1
		
		if np.max(np.abs(error)) < tol:
		
			xrpo[gindices, :] = xi
			tofrpo[gindices] = Ti
			phirpo[gindices] = phii
			converged = True
		
		if iteration == itermax or AdaptiveFail:

			if AdaptiveFail:
				
				print "Adaptive step failed, looks kinda hopeless."
			
			print "did not converged in given maximum number of steps"
			print "exitting..."
			xrpo[gindices, :] = xi
			tofrpo[gindices] = Ti
			phirpo[gindices] = phii
			break
		
	
	print "x_rpo:"
	print xrpo
	print "T_rpo:"
	print tofrpo
	print "phi_rpo:"
	print phirpo
	
	plot(np.arange(iteration), earray)
	xlabel('Iteration')
	ylabel('Error')
	savefig('image/errorrposearch'+str(i)+'.png', bbox_inches='tight', dpi=150)
	#show()
	
	#raw_input("Press Enter to continue...")

np.savetxt('data/tofrpo.dat',tofrpo)
np.savetxt('data/xrpo.dat', xrpo)
np.savetxt('data/phirpo.dat', phirpo)

periods = np.zeros(np.max(group))

for i in range(1,int(np.max(group))+1):
		#Get corresponding indices for the group i:
		indices = np.argwhere(group == i)
		#Reshape them in a 1D array:
		indices = indices.reshape(np.size(indices))
		periods[i-1] = np.sum(tofrpo[indices])


np.savetxt('data/periods.dat', periods)

sys.exit("Bye!")
