#
# reconstructfull.py
#

from __future__ import unicode_literals
import KS
import numpy as np
from scipy.linalg import expm

tauxhatphit = np.loadtxt("data/solutionscaledtime.dat") #Load the scaled time solution

xhattau = tauxhatphit[:,1:np.size(tauxhatphit,1)-2]
phitau = tauxhatphit[:,np.size(tauxhatphit,1)-2]
ttau = tauxhatphit[:,np.size(tauxhatphit,1)-1]

Nm1 = int((np.size(xhattau,1))/2+1) - 1 	# N of modes -1
T = KS.generator(Nm1) # Lie algebra generator

def LieElement(phi):
	g = expm(phi*T)
	return g

skip = 1
xtau = np.array([np.dot(LieElement(phitau[i]), xhattau[i,:]) for i in range(0, np.size(phitau,0), skip)], float)
#print xtau

tx = np.array([ttau[range(0, np.size(phitau,0), skip)]], float).transpose() #Add time to the first collumn of tx
tx = np.append(tx, xtau, axis=1) #Add solution to the following collumns of tx

np.savetxt("data/solutionreconstructed.dat", tx)

# Plot the solution:
#Import plotting modules:
from pylab import figure, grid, hold, legend, savefig
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

fig = plt.figure() #Create a figure instance
ax = fig.gca(projection = "3d") #Create an axis instance

ax.plot(xtau[:,0], xtau[:,1], xtau[:,2], lw=0.5)
ax.set_xlabel('$x_1$', fontsize=18)
ax.set_ylabel('$y_1$', fontsize=18)
ax.set_zlabel('$x_2$', fontsize=18)
ax.view_init(20,120)
savefig("image/reconstructedsolution.png", bbox_inches='tight', dpi=100)

plt.tight_layout()
plt.show()
