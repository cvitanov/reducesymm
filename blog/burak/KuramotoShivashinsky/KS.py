#
# twomode.py
#
"""
This module holds the definitions of velocity functions, Lie algebra 
generators and Lie algebra elements specific to the Two mode system
"""

import numpy as np

def vfullssp(x, t):
    """
    Velocity function in the full state space.
	
    Arguments:
    x : vector of mode values
		v = [Re[a1] Im[a1] Re[a2] Im[a2] ... Re[aN-1], Im[aN-1]]
    N : number of modes to be used (0 to N-1)
    """
	
    L = 22
    Ltilde = 22/np.pi
    N = int(np.size(x)/2+1)
    #From state space variables to complex modes:
    a = np.zeros(N, complex)
    a[1:N] = np.array([x[2*i] + 1j*x[(2*i) + 1] for i  in range(N-1)], complex)
    #Velocity function:
    adot = np.zeros(N,complex)
    #adot[0] = 0, range starts from k = 1
    for k in range(1,N):
	
	qk = complex(k)/Ltilde
	
	nonlinearterm = complex(0 + 1j*0)
    	
	for m in range(-N+1, N):
	    if m < 0:
		am = np.conjugate(a[-m])
	    else:
		am = a[m]
	    
	    if (k-m) in range(-N+1, 0):
		akMinm = np.conjugate(a[m-k])		
	    elif (k-m) in range(0, N):
		akMinm = a[k-m]	    
	    else:
		akMinm = 0
	    
	    nonlinearterm += am*akMinm
	 
	adot[k] = (qk**2-qk**4)*a[k] - 1j*(qk/2)*nonlinearterm    
    
    vel = np.zeros((N-1)*2, float)
    vel[range(0,(N-1)*2-1, 2)] = np.array(np.real(adot[1:N]) ,float)
    vel[range(1,(N-1)*2, 2)] = np.array(np.imag(adot[1:N]) ,float)
    
    return vel

def generator(Nm1):
	"""
	Generator of infinitesimal SO(2) transformations
	Nm1 : number of modes counted from 1st mode
	"""
	
	T = np.zeros((2*Nm1,2*Nm1), float)
    
	for i in range(Nm1):
		T[2*i, 2*i+1] = i+1
		T[2*i+1, 2*i] = -i-1
		
	return T

def vscaledtime(x, tau):
	"""
	Velocity function for the symmetry reduced two mode system where  time 
	is scaled as 
	dt = x1 dtau
	everything is within the 1st mode slice of:
	xhat = (1,0,0,0)
	"""
	
	Nm1 = int(np.size(x)/2+1) - 1 
	
	T = generator(Nm1)
	vhat = x[0]*np.array(vfullssp(x, tau)) + vfullssp(x, tau)[1]*np.dot(T,x)
	#print vhat 
	return vhat

def vphasescaledtime(x, tau):
    
    vphase=-vfullssp(x, tau)[1]
    
    return vphase

def dtdtau(x, tau):
    
    dtdtauscaledtime = x[0]
    
    return dtdtauscaledtime

def LieElement(N,phi):
    """
    Lie Element of SO(2) transformations
    N : number of modes counted from 1st mode
    """
    
    T = generator(N)
    g = expm(phi*T)

    return g


def StabilityMatrix(N,x):
	
	x1, x2, y1, y2 = x
	
	A = np.zeros((2*N,2*N), float)
		
	return A
