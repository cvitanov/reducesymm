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

Zeta0 = 1
Trace00 = 0

s, z = sympy.symbols('s,z')

for rpono in range(1,26):
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
        
    Zeta0 = Zeta0 * (1 - (sympy.exp(-s*T)*z**TopLength)/np.abs(floquet[0]))

conn.close()

from sympy import degree, LT
#Zeta0 = sympy.expand(Zeta0)    
#Zeta00 = []
#ExpOrder = 4

raw_input("fdasfdsafdsa")
for j in range(ExpOrder,0,-1):
    while degree(Zeta00) > j:
    
        Zeta0 = Zeta0 - LT(Zeta0)
        if degree(Zeta0) <= ExpOrder:
            Zeta00.insert(0, (np.float(Zeta00.subs(z,1))))


#print Zeta00
#print Zeta00.subs(z, 1)
print Zeta00Val
plot(range(1, ExpOrder+1),Zeta00Val)
plt.show()
