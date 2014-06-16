import numpy as np
import sqlite3
from StringIO import StringIO 
from subprocess import call
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
Ncycle = 26

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

z = sympy.symbols('z')
Trace = 0
Nexpansion = 9
#Compute trace according to ChaosBook v14 eq. 20.13
for i in range(Ncycle):
    ni = rpos[i][0]
    Ti = rpos[i][1]
    Lambda = rpos[i][2]
    Sum = 0
    r = 1
    while ni*r <= Nexpansion:
        Sum = Sum + (z**(ni*r))/np.abs(1-Lambda**r)
        r += 1
    
    Trace = Trace + Ti * Sum
#Extract coefficients of the expansion:
C = []
for i in range(Nexpansion):
    C.insert(0, sympy.LC(Trace))
    Trace = Trace - sympy.LT(Trace)
#Compute Coefficients of the determinant expansion according to ChaosBook v14 eq. 20.15    
Q = [C[0]]
SpectDetPoly = np.array([-Q[0], 1], float)
F00 = [1 - sum(Q)]
for n in range(2, Nexpansion+1):
    Qn = C[n-1]
    for k in range(1,n):
        #print n,k
        Qn = Qn - C[n - k - 1]*Q[k - 1] 
    Qn = Qn / float(n)
    Q.append(Qn)
    F00.append(1 - sum(Q))    
    SpectDetPoly = np.insert(SpectDetPoly, 0, -Qn)
