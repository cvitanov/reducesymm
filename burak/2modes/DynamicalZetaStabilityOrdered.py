import numpy as np
import sqlite3
from StringIO import StringIO 
from subprocess import call
import sys
import sympy
from scipy.optimize import fsolve
#Initiate plotting environment:
import matplotlib as mpl
from pylab import plot, xlabel, ylabel, show, savefig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#ExpOrder = 8
EscapeRate = []
ConservationRule = []
Ncycle = 22
Lambdath = 1e-6 #Stability Threshold

for ExpOrder in range(1,8):

    conn = sqlite3.connect('data/rpoall.db')
    c = conn.cursor()
    
    Zeta0 = 1
    s, z = sympy.symbols('s,z')
    
    for rpono in range(1,Ncycle+1):
        c.execute("SELECT * FROM rpos WHERE rpono = "+str(rpono))
        a = c.fetchall()
        
        rpono,  itinerary, x1, y1, x2, y2, T, phifetch , floquetFetch= a[0]
        #Format the Floquet exponents such that it would be readable by loadtxt as
        #a numpy array
        floquetobj = StringIO(str(floquetFetch).strip('[ ]'))
        floquet = np.loadtxt(floquetobj)
        TopLength = len(str(itinerary))
        #if str(itinerary)=='001':
            #continue
            
        Zeta0 = Zeta0 * (1 - (sympy.exp(-s*T)/np.abs(floquet[0]))
        Zeta0 = Zeta0.expand()
        while sympy.degree(Zeta0, z) > ExpOrder:
            Zeta0 = Zeta0 - sympy.LT(Zeta0, z)
            Zeta0 = sympy.collect(Zeta0, z)
    
    conn.close()
    
    Zeta0 = Zeta0.subs(z,1)
    
    ConservationRule.append(Zeta0.subs(s,0))
    print "Conservation rule:", ConservationRule
    
    f = sympy.lambdify(s, Zeta0, "numpy")
    def fcomplex(xy):
        x, y = xy
        fcomp = np.array([np.real(f(x + 1j*y)), np.imag(f(x + 1j*y))])
        return fcomp
    #found = False
    #step = 1e-4
    #splusnext = 0
    #sminusnext = 0
    #while not(found):
        #splus = splusnext
        #sminus = sminusnext
        #splusnext += step
        #sminusnext -= step
        #if Zeta0.subs(s, sminus)*Zeta0.subs(s, sminusnext) < 0:
            #EscapeRate.append(float(sympy.nsolve(Zeta0, (sminus+sminusnext)/2.0, tol=1e-7)))
            #found = True
        #if Zeta0.subs(s, splus)*Zeta0.subs(s, splusnext) < 0:
            #EscapeRate.append(float(sympy.nsolve(Zeta0, (splus+splusnext)/2.0, tol=1e-7)))
            #found = True
    #EscapeRate.append(fsolve(fcomplex, [0, 0]))        
    EscapeRate.append(-float(fsolve(f, 0)))
    print "EscapeRate:", EscapeRate
    
escrate = np.array(EscapeRate, float)
np.savetxt('data/escratedynzeta.dat', escrate)
