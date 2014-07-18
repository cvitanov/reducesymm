import numpy as np
import sqlite3
from StringIO import StringIO 
from subprocess import call
from scipy.optimize import fsolve
import sys
#import sympy
#Initiate plotting environment:
import matplotlib as mpl
from pylab import plot, xlabel, ylabel, show, savefig
import matplotlib.pyplot as plt
import twomode

FiniteGrammar=True

if FiniteGrammar:
    Ncycle = 79
    NmaxExp = 12
    conn = sqlite3.connect('data/rpo.db')

else:
    Ncycle = 133
    NmaxExp = 12
    conn = sqlite3.connect('data/rpoall.db')

c = conn.cursor()
rpos = []

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
    rpos.append([TopLength, T, floquet[0], phi])

conn.close()

orig_stdout = sys.stdout

f = file('data/rpotext.dat', 'w')
sys.stdout = f

for i in range(len(rpos)):
    print str(rpos[i]).strip("[]").replace(",", "")

sys.stdout = orig_stdout
f.close()
