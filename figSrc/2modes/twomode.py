#
# twomode.py
#
"""
This module holds the definitions of velocity functions, Lie algebra 
generators and Lie algebra elements specific to the Two mode system
"""

import numpy as np

def vinvpol(x, t, p):
    """
    Velocity function for the invariant polynomial basis.
    
    Arguments:
    x : vector of invariant polynomials
    	v = [u,v,w,q]
    t : time
    p : vector of parameters:
    	p = [mu1, a1, b1, c1, mu2, a2, b2, c2, e2]
    """
	
    u, v, w, q = x
    mu1, a1, b1, c1, e1, mu2, a2, b2, c2, e2 = p

    #The velocity function v = d(u,v,w,q)/dt:
    vel = [2*mu1*u + 2*a1*u**2 + 2*b1*u*v + c1*w,
    	 2*mu2*v + 2*a2*u*v + 2*b2*v**2 + c2*w,
    	(2*mu1+mu2)*w + (2*a1+a2)*u*w + (2*b1+b2)*v*w + 4*c1*u*v+2*c2*u**2-e2*q,
    	(2*mu1 + mu2)*q + (2*a1 + a2)*u*q + (2*b1 + b2)*v*q + e2*w]

    return vel

def vfullssp(x, t, p):
    """
    Velocity function in the full state space.
	
    Arguments:
    x : vector of invariant polynomials
    	v = [u,v,w,q]
    t : time
    p : vector of parameters:
    	p = [mu1, a1, b1, c1, mu2, a2, b2, c2, e2]
    """
	
    #x1, x2, y1, y2 = x
    x1, y1, x2, y2 = x
    r1sq = x1**2 + y1**2
    r2sq = x2**2 + y2**2
	
    #mu1, a1, b1, c1, mu2, a2, b2, c2, e2 = p
    mu1, a1, b1, c1, e1, mu2, a2, b2, c2, e2 = p
    
    #The velocity function v = d(x1,x2,y1,y2)/dt:
    #vel = [mu1*x1 + a1*x1**3 + b1*x1*y1**2  + c1*x1*y1 + a1*x1*x2**2 + b1*x1*y2**2 + c1*x2*y2,
    #	   mu1*x2 + a1*x1**2*x2 + c1*x1*y2 + b1*y1**2*x2 - c1*y1*x2 + a1*x2**3 + b1*x2*y2**2,
    #	   mu2*y1 + a2*x1**2*y1 + c2*x1**2 + b2*y1**3 + a2*y1*x2**2 + b2*y1*y2**2 - c2*x2**2 + e2*y2,
    #	   mu2*y2 + a2*x1**2*y2 + 2*c2*x1*x2 + b2*y1**2*y2 - e2*y1 + a2*x2**2*y2 + b2*y2**3]
    vel = [(mu1 + c1*x2 + a1*r1sq + b1*r2sq)*x1 + c1*y1*y2 + e1*y1,
    	   (mu1 - c1*x2 + a1*r1sq + b1*r2sq)*y1 + c1*y2*x1 - e1*x1,
    	   (mu2 + a2*r1sq + b2*r2sq)*x2 + c2*(x1**2 - y1**2) + e2*y2,
    	   (mu2 + a2*r1sq + b2*r2sq)*y2 + 2*c2*x1*y1 - e2*x2]

    return vel

def generator():
    """
    Generator of infinitesimal SO(2) transformations for the two mode system
    """
    T = np.array([[0,-1,0,0],
    		  [1,0,0,0],
    		  [0,0,0,-2],
    		  [0,0,2,0]], 
    			 float)

    return T

def LieElement(phi):
    """
    Lie Element of SO(2) transformations for the two mode system
    """
    g = np.array([[np.cos(phi),-np.sin(phi),0,0],
    			 [np.sin(phi),np.cos(phi),0,0],
    			 [0,0,np.cos(2*phi),-np.sin(2*phi)],
    			 [0,0,np.sin(2*phi),np.cos(2*phi)]],
    			 float)

    return g

def vscaledtime(x, tau, p):
	"""
	Velocity function for the symmetry reduced two mode system where  time 
	is scaled as 
	dt = x1 dtau
	everything is within the 1st mode slice of:
	xhat = (1,0,0,0)
	"""
	
	T = generator()
	vhat = x[0]*np.array(vfullssp(x, tau, p)) - vfullssp(x, tau, p)[1]*np.dot(T,x)
	#print vhat 
	return vhat
	
def StabilityMatrix(x, p):
	
	#x1, x2, y1, y2 = x
	x1, y1, x2, y2 = x
	mu1, a1, b1, c1, e1, mu2, a2, b2, c2, e2 = p
	
	#A = np.array([[3*x1**2*a1 + x2**2*a1 + y1**2*b1 + y2**2*b1 + y1*c1 + mu1, 2*x1*x2*a1 + y2*c1, 2*x1*y1*b1 + x1*c1, 2*x1*y2*b1 + x2*c1],
				  #[2*x1*x2*a1 + y2*c1, x1**2*a1 + 3*x2**2*a1 + y1**2*b1 + y2**2*b1 - y1*c1 + mu1, 2*x2*y1*b1 - x2*c1, 2*x2*y2*b1 + x1*c1],
				  #[2*x1*y1*a2 + 2*x1*c2, 2*x2*y1*a2 - 2*x2*c2, x1**2*a2 + x2**2*a2 + 3*y1**2*b2 + y2**2*b2 + mu2, 2*y1*y2*b2+e2],
				  #[2*x1*y2*a2 + 2*x2*c2, 2*x2*y2*a2 + 2*x1*c2, 2*y1*y2*b2 - e2, x1**2*a2 + x2**2*a2 + y1**2*b2 + 3*y2**2*b2  + mu2]])
				  
	A = np.array([[3*x1**2*a1 + y1**2*a1 + x2**2*b1 + y2**2*b1 + x2*c1 + mu1, 2*x1*y1*a1 + y2*c1, 2*x1*x2*b1 + x1*c1, 2*x1*y2*b1 + y1*c1],
				  [2*x1*y1*a1 + y2*c1, x1**2*a1 + 3*y1**2*a1 + x2**2*b1 + y2**2*b1 - x2*c1 + mu1, 2*y1*x2*b1 - y1*c1, 2*y1*y2*b1 + x1*c1],
				  [2*x1*x2*a2 + 2*x1*c2, 2*y1*x2*a2 - 2*y1*c2, x1**2*a2 + y1**2*a2 + 3*x2**2*b2 + y2**2*b2 + mu2, 2*x2*y2*b2+e2],
				  [2*x1*y2*a2 + 2*y1*c2, 2*y1*y2*a2 + 2*x1*c2, 2*x2*y2*b2 - e2, x1**2*a2 + y1**2*a2 + x2**2*b2 + 3*y2**2*b2  + mu2]])
	
	return A

def velocityvar(x, t, p):
	"""
	Velocity function for the variational equation.
	Takes 4+4x4-D input x and parameter vector p as the input and gives
	4+4x4-D velocity output
	Relevant notes: %http://www.cs.colorado.edu/~lizb/chaos/variational-notes.pdf
	"""
	
	vel = np.zeros(4+4*4) #Generate the dummy vector
	vel[0:4] = vfullssp(x[0:4],t,p) #First 4 elements are the regular velocity
	mvars = x[4:4+4*4].reshape(4,4) #Read the variations matrix
	mvelvar = np.dot(StabilityMatrix(x[0:4],p), mvars) #Velocity of variations in matrix form
	
	vel[4:4+4*4] = mvelvar.reshape(16) 
	
	return vel

def ssp2invpol(x):
	"""
	Takes statespace coordinates and returns corresponding point in the
	invariant polynomial bases
	"""
	
	x1, y1, x2, y2 = x
	#Convert to complex modes:
	z1 = x1 + 1j*y1 
	z2 = x2 + 1j*y2
	#Convert to inv pols:
	u = z1*np.conj(z1)
	v = z2*np.conj(z2)
	w = (z1**2)*np.conj(z2) + (np.conj(z1)**2)*z2
	q = ((z1**2)*np.conj(z2) - (np.conj(z1)**2)*z2)/1j
	
	uu = np.array(np.real([u, v, w, q]), float)
	return uu

def ssp2sspRed2(xsol):
	"""
	Takes the solution in the full statespace and outputs a symmetry reduced
	representation where the 2nd mode phase is fixed.
	"""
	#Construct complex solution vector
	asol = np.array([xsol[:,0] + 1j*xsol[:,1], xsol[:,2] + 1j*xsol[:,3]], 
	complex).transpose() 
	#Construct phi_2(t)
	phi = np.array([- np.angle(xsol[i, 2] + xsol[i, 3]*1j)/2 for i in 
	range(np.size(xsol, 0))], float)
	#Construct g(phi(t))
	gUslice = np.array([[np.exp( 1j*phi[i]), np.exp( 1j*2*phi[i])]  for 
	i in range(np.size(xsol, 0))], complex)
	asolhat = np.multiply(gUslice, asol)
	#Phase doubling:
	ephitilde = np.array([[np.exp( 1j*np.angle(asolhat[i,0])), np.exp(0)]  
	for i in range(np.size(xsol, 0))], complex)
	asoltilde = np.multiply(ephitilde, asolhat)
	#Back to real rep:
	xsoltilde = np.array([np.real(asoltilde[:,0]), np.imag(asoltilde[:,0]), 
	np.real(asoltilde[:,1]), np.imag(asoltilde[:,1])], float).transpose()
	return xsoltilde
