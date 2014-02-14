#!/usr/bin/env python
"""
http://www.cs.yale.edu/publications/techreports/tr1435.pdf
A Highly Accurate Solver for Stiﬀ Ordinary Diﬀerential
Equations
Dan Kushnir∗ Vladimir Rokhlin†
"""
import numpy as np

def KushnirRokhlin(f, x0, t)
    """
    USAGE:
        x = KushnirRokhlin(f, x0, t)

    INPUT:
        f     - function of x and t equal to dx/dt.  x may be multivalued,
                in which case it should a list or a NumPy array.  In this
                case f must return a NumPy array with the same dimension
                as x.
        x0    - the initial condition(s).  Specifies the value of x when
                t = t[0].  Can be either a scalar or a list or NumPy array
                if a system of equations is being solved.
        t     - list or NumPy array of t values to compute solution at.
                t[0] is the the initial condition point, and the difference
                h=t[i+1]-t[i] determines the step size h.

    OUTPUT:
        x     - NumPy array containing solution values corresponding to each
                entry in t array.  If a system is being solved, x will be
                an array of arrays.
    """

    n = len( t )
    x = numpy.array( [x0] * n )
    for i in xrange( n - 1 ):
        x[i+1] = x[i] + ( t[i+1] - t[i] ) * f( x[i], t[i] )

    return x
