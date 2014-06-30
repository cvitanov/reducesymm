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
import matplotlib.pyplot as plt
import twomode

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
    
    c.execute("SELECT phi FROM rpo"+str(itinerary))
    a = c.fetchall()
    phi = -sum([float(str(a[i]).strip('( ,)')) for i in range(len(a))])
    
    tzeta = 1.0/abs(floquet[0])
    
    rpos.append([TopLength, T, floquet[0], phi, tzeta])

conn.close()

