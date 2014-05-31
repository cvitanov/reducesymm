import numpy as np
import sqlite3
#Initiate plotting environment:
import matplotlib as mpl
from pylab import plot, xlabel, ylabel, show, savefig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from subprocess import call
#Import twomode system module
import twomode 

conn = sqlite3.connect('data/rpo.db')
c = conn.cursor()

for rpono in range(1,3):
    c.execute("SELECT * FROM rpos WHERE rpono = "+str(rpono))
    a = c.fetchall()
    
    rpono,  itinerary, x1, y1, x2, y2, T, phi = a[0]
    z1 = x1 + 1j*y1
    z2 = x2 + 1j*y2
    slicePhase = np.angle(z1)
    z1 = np.exp(-1j*slicePhase)*z1
    z2 = np.exp(-2*1j*slicePhase)*z2
    
    rpo = [np.real(z1), np.imag(z1), np.real(z2), np.imag(z2), 0]
    dt = 0.01
    tsol = np.linspace(0, T, np.floor(T/dt)+1)
    xphisol = twomode.intslice(rpo, tsol)
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #Modify axis colors:
    ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
    ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
    ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
            
    ax.plot(xphisol[:,0], 
    xphisol[:,2], 
    xphisol[:,3], linewidth=1, color='#3c5f96')
    
    ax.set_xlabel('\n $\hat{x}_1$ \t  ', fontsize=32)
    ax.set_ylabel('\n $\hat{x}_2$ \t', fontsize=32)
    ax.set_zlabel('$\hat{y}_2$   ', fontsize=32)
    
    savefig('image/'+str(itinerary)+'.png', bbox_inches='tight', dpi=100) 
    
conn.close()    
