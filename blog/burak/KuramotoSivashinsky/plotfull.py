#
# plotfull.py
#

from __future__ import unicode_literals
import KS
import numpy as np
from scipy.linalg import expm

tx = np.loadtxt("data/sspsolution.dat") #Load the scaled time solution

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

ax.plot(tx[:,1], tx[:,2], tx[:,3], lw=0.5)
ax.set_xlabel('$x_1$', fontsize=18)
ax.set_ylabel('$y_1$', fontsize=18)
ax.set_zlabel('$x_2$', fontsize=18)
ax.view_init(20,120)
savefig("image/ssp.png", bbox_inches='tight', dpi=100)

plt.tight_layout()
plt.show()
