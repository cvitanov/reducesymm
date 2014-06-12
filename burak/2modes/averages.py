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

escapeRate = 1

z = sympy.symbols('z')

for rpono in range(1,19):
    c.execute("SELECT * FROM rpos WHERE rpono = "+str(rpono))
    a = c.fetchall()
    
    rpono,  itinerary, x1, y1, x2, y2, T, phifetch , floquetFetch= a[0]
    #Format the Floquet exponents such that it would be readable by loadtxt as
    #a numpy array
    floquetobj = StringIO(str(floquetFetch).strip('[ ]'))
    floquet = np.loadtxt(floquetobj)
    TopLength = len(str(itinerary))
    #if str(itinerary)=='001':
    #    continue
        
    escapeRate = escapeRate * (1 - (z**TopLength)/np.abs(floquet[0]))

from sympy import degree, LT
escapeRate = sympy.expand(escapeRate)    
escapeRateVal = []
ExpOrder = 8
for j in range(ExpOrder,0,-1):
    while degree(escapeRate) > j:
    
        escapeRate = escapeRate - LT(escapeRate)
        if degree(escapeRate) <= ExpOrder:
            escapeRateVal.insert(0, np.abs(escapeRate.subs(z,1)))


#print escapeRate
#print escapeRate.subs(z, 1)
print escapeRateVal
plot(range(1, ExpOrder+1),escapeRateVal)
plt.show()
