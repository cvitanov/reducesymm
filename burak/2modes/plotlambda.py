import numpy as np
import sqlite3
#Initiate plotting environment:
import matplotlib as mpl
from pylab import plot, xlabel, ylabel, show, savefig, xlim, ylim
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from subprocess import call
import sys
#Import twomode system module
import twomode 
from StringIO import StringIO 
conn = sqlite3.connect('data/rpoall.db')
c = conn.cursor()

Ncycle = 127

for rpono in range(1,Ncycle+1):
    c.execute("SELECT * FROM rpos WHERE rpono = "+str(rpono))
    a = c.fetchall()
    
    rpono,  itinerary, x1, y1, x2, y2, T, phifetch , floquetFetch= a[0]
    #Format the Floquet exponents such that it would be readable by loadtxt as
    #a numpy array
    floquetobj = StringIO(str(floquetFetch).strip('[ ]'))
    floquet = np.loadtxt(floquetobj)
    
    TopLength = len(str(itinerary))
    
    fig = plt.figure(2, figsize=(8,6))
    plot(np.log(np.abs(floquet[0]))/T ,1.0/TopLength, '.', ms=16)
    
    plt.hold(True)
    plt.grid(True)

ylim(0,0.6)
ax = fig.gca()
Nticks = 8
xticks = np.linspace(0.05, 0.40, Nticks)
ax.set_xticks(xticks) 
ax.set_xticklabels(["$%.2f$" % xtik for xtik in xticks], fontsize=12); 

yticks = np.array([1.0/n for n in range(2,13)], float)
ax.set_yticks(yticks)
ax.set_yticklabels(["%s" % n for n in ['$1/2$', '$1/3$', '$1/4$', '', '$1/6$', '', '$1/8$', 
                                        '', '$1/10$', '', '\n $1/12$']], fontsize=14); 
xlabel('$\lambda$', fontsize=24)
ylabel('$1/n$', fontsize=24)

savefig('lambdaDist.pdf', bbox_inches='tight', dpi=100) 
conn.close()    

plt.show()
