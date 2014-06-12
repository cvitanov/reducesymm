import numpy as np
import sqlite3
#Initiate plotting environment:
import matplotlib as mpl
from pylab import plot, xlabel, ylabel, show, savefig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from subprocess import call
import sys
#Import twomode system module
import twomode 
from StringIO import StringIO 
conn = sqlite3.connect('data/rpo.db')
c = conn.cursor()

f = open("tex/twomoderpos.tex", "w")

f.write("\\begin{table}\n")
f.write("\t\\begin{tabular}{c|c|c|c|c}\n")

f.write("\tItinerary & $(x_{1,RPO}, y_{1,RPO}, x_{2,RPO}, y_{2,RPO})$ & Period & Phase Shift & Floquet Multipliers \\\\ \n")
f.write("\t\\hline\n")

for rpono in range(1,27):
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
    phi = xphisol[-1,-1]
    
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
    
    f.write("\t%s & " % itinerary)
    f.write("(%5.8f, " % float(rpo[0])) 
    f.write("%5.1f, " % abs(float(rpo[1]))) 
    f.write("%5.8f, " % float(rpo[2])) 
    f.write("%5.8f) & " % float(rpo[3])) 
    f.write("%5.8f & " % float(T))
    f.write("%5.8f & " % float(phi))
    f.write("(%5.8f, " % float(floquet[0])) 
    f.write("%5.8f, " % float(floquet[1]))
    f.write("%5.8f, " % float(floquet[2])) 
    f.write("%0.4g) \\\\ \n " % float(floquet[3])) 

f.write("\t\\end{tabular}\n")
f.write("\t\\caption{\\rpo s of the \\twoMode\\ system. \
Parameter values \\reftab{tab:pars}\,(a).}\n")
f.write("\t\\label{t-rpo2modeupto8}\n")
f.write("\\end{table}")

f.close()

#plt.show()
conn.close()    
