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
    mu1, a1, b1, c1, mu2, a2, b2, c2, e2 = p

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
	
    x1, x2, y1, y2 = x
    mu1, a1, b1, c1, mu2, a2, b2, c2, e2 = p

    #The velocity function v = d(u,v,w,q)/dt:
    vel = [mu1*x1 + a1*x1**3 + b1*x1*y1**2  + c1*x1*y1 + a1*x1*x2**2 + b1*x1*y2**2 + c1*x2*y2,
    	   mu1*x2 + a1*x1**2*x2 + c1*x1*y2 + b1*y1**2*x2 - c1*y1*x2 + a1*x2**3 + b1*x2*y2**2,
    	   mu2*y1 + a2*x1**2*y1 + c2*x1**2 + b2*y1**3 + a2*y1*x2**2 + b2*y1*y2**2 - c2*x2**2 + e2*y2,
    	   mu2*y2 + a2*x1**2*y2 + 2*c2*x1*x2 + b2*y1**2*y2 - e2*y1 + a2*x2**2*y2 + b2*y2**3]

    return vel

def generator():
    """
    Generator of infinitesimal SO(2) transformations for the two mode system
    """
    T = np.array([[0,1,0,0],
    			 [-1,0,0,0],
    			 [0,0,0,2],
    			 [0,0,-2,0]], 
    			 float)

    return T

def LieElement(phi):
    """
    Lie Element of SO(2) transformations for the two mode system
    """
    g = np.array([[np.cos(phi),np.sin(phi),0,0],
    			 [-np.sin(phi),np.cos(phi),0,0],
    			 [0,0,np.cos(2*phi),np.sin(2*phi)],
    			 [0,0,-np.sin(2*phi),np.cos(2*phi)]],
    			 float)

    return g
