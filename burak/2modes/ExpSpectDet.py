import numpy as np
import sqlite3
from StringIO import StringIO 
from subprocess import call
from scipy.optimize import fsolve
import sys
import sympy
#Initiate plotting environment:
import matplotlib as mpl
from pylab import plot, xlabel, ylabel, show, savefig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

conn = sqlite3.connect('data/rpo.db')
c = conn.cursor()
rpos = []
Ncycle = 36

for rpono in range(1,Ncycle + 1):
    c.execute("SELECT * FROM rpos WHERE rpono = "+str(rpono))
    a = c.fetchall()
    
    rpono,  itinerary, x1, y1, x2, y2, T, phifetch , floquetFetch= a[0]
    #Format the Floquet exponents such that it would be readable by loadtxt as
    #a numpy array
    floquetobj = StringIO(str(floquetFetch).strip('[ ]'))
    floquet = np.loadtxt(floquetobj)
    TopLength = len(str(itinerary))
    rpos.append([TopLength, T, floquet[0]])

conn.close()

s, z = sympy.symbols('s z')
EscapeRate = []
ConservationRule = []

for Nexpansion in range(1,11):
    SpectralDeterminant = 1
    #Nexpansion = 3
    Exponent = 0
    #Compute Spectral Determinant according to ChaosBook v14 eq. 19.6
    for i in range(Ncycle):
        npr = rpos[i][0]
        Tp = rpos[i][1]
        Lambda = rpos[i][2]
        Sum = 0
        r = 1
        while npr*r <= Nexpansion:
            Sum = Sum + (sympy.exp(-r*s*Tp)*z**(npr*r))/(r*abs(1.0-Lambda**r))
            r += 1
        
        ExpSum = sympy.series(sympy.exp(-Sum), z, n=Nexpansion+1).subs(sympy.O(z**(Nexpansion+1)), 0)
        SpectralDeterminant = (SpectralDeterminant * ExpSum).expand()
        while sympy.degree(SpectralDeterminant, z) > Nexpansion:
            SpectralDeterminant = SpectralDeterminant - sympy.LT(SpectralDeterminant, z)
            SpectralDeterminant = sympy.collect(SpectralDeterminant, z)
    
    SpectralDeterminant = SpectralDeterminant.subs(z,1)    
        #Exponent = Exponent - Sum
    #Extract coefficients of the expansion:
    #SpectralDeterminant = sympy.series(sympy.exp(Exponent),z,n=Nexpansion+1).subs(sympy.O(z**(Nexpansion+1)), 0)
    #SpectralDeterminant = sympy.series(sympy.exp(Exponent),z,n=Nexpansion+1).subs(sympy.O(z**(Nexpansion+1)), 0).subs(z, 1)
    ConservationRule.append(SpectralDeterminant.subs(s,0))
    #EscapeRate.append(float(sympy.nsolve(SpectralDeterminant, 0)))
    f = sympy.lambdify(s, SpectralDeterminant, "numpy")
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
        #if f(sminus)*f(sminusnext) < 0:
            #EscapeRate.append(float(fsolve(f, (sminus + sminusnext)/2.0)))
            #found = True
        #if f(splus)*f(splusnext) < 0:
            #EscapeRate.append(float(fsolve(f, (splus + splusnext)/2.0)))
            #found = True
    EscapeRate.append(float(fsolve(f, 0)))
    #EscapeRate.append(fsolve(fcomplex, [0, 0]))
    print EscapeRate


f = open("tex/s0.tex", "w")

f.write("\\begin{table}\n")
f.write("\t\\begin{tabular}{c|c}\n")

f.write("\t N & s_0 \\\\ \n")
f.write("\t\\hline\n")

for i in range(len(EscapeRate)):
    f.write("\t%s & " % str(i+1))
    f.write("%5.8f \\\\ \n " % float(EscapeRate[i]))
    
f.write("\t\\end{tabular}\n")
f.write("\t\\caption{Leading zero of the spectral determinant (\\beta = 0) \
computed using the finite grammar approximation}\n")
f.write("\t\\label{t-s0}\n")
f.write("\\end{table}")

f.close()     

escrate = np.array(EscapeRate, float)
np.savetxt('data/escratespecdet.dat', escrate)
