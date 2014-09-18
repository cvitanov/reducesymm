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

Ncycle = 4

fig = plt.figure()
ax = fig.gca(projection='3d')
#Modify axis colors:
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))

plt.hold(True)
    
for rpono in range(1,Ncycle+1):
    c.execute("SELECT * FROM rpos WHERE rpono = "+str(rpono))
    a = c.fetchall()
    
    rpono,  itinerary, x1, y1, x2, y2, T, phifetch , floquetFetch= a[0]
    #Format the Floquet exponents such that it would be readable by loadtxt as
    #a numpy array
    floquetobj = StringIO(str(floquetFetch).strip('[ ]'))
    floquet = np.loadtxt(floquetobj)
    
    z1 = x1 + 1j*y1
    z2 = x2 + 1j*y2
    slicePhase = np.angle(z1)
    z1 = np.exp(-1j*slicePhase)*z1
    z2 = np.exp(-2*1j*slicePhase)*z2
    
    rpo = [np.real(z1), np.imag(z1), np.real(z2), np.imag(z2), 0]
    dt = 0.01
    tsol = np.linspace(0, T, np.floor(T/dt)+1)
    xphisol = twomode.intslice(rpo, tsol)
    #print phifetch, xphisol[-1,-1]
    #phi = xphisol[-1,-1]
    phi = -phifetch
    
    TopLength = len(str(itinerary))
    
    if str(itinerary) == '001':
        lw = 3
    else:
        lw=2
            
    ax.plot(xphisol[:,0], 
    xphisol[:,2], 
    xphisol[:,3], linewidth=lw) #, color='#3c5f96')
    
    ax.set_xlabel('\n $\hat{x}_1$ \t  ', fontsize=32)
    ax.set_ylabel('\n $\hat{x}_2$ \t', fontsize=32)
    ax.set_zlabel('$\hat{y}_2$   ', fontsize=32)
    
conn.close()    

ax.view_init(15,30)

Nticks = 3

ax.set_xlabel('\n $\hat{x}_1$ \t  ', fontsize=32)
ax.set_ylabel('\n $\hat{x}_2$ \t', fontsize=32)
ax.set_zlabel('$\hat{y}_2$   ', fontsize=32)

xticks = np.linspace(0, 2, Nticks)
ax.set_xticks(xticks) 
ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=24); 

yticks = np.linspace(-1.8, 0, Nticks)
ax.set_yticks(yticks)
ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=24); 

zticks = np.linspace(-0.4, 0.4, Nticks)
ax.set_zticks(zticks)
#ax.set_zticklabels(["$%.1f$ " % ztik for ztik in zticks], fontsize=24); 
ax.set_zticklabels(["$-0.4$ \t ", " $0.0$  ", " $0.4$  "], fontsize=24); 

savefig('rpofirst'+str(Ncycle)+'.pdf', bbox_inches='tight', dpi=100) 

plt.show()

