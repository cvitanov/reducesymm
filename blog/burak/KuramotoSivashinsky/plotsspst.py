#
# plotfull.py
#

from __future__ import unicode_literals
import KS
import numpy as np
from scipy.linalg import expm

tx = np.loadtxt("data/solutionreconstructed.dat") #Load the scaled time solution

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
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
#ax.view_init(20,-65)
ax.view_init(25,120)
savefig("image/sspst.png", bbox_inches='tight', dpi=100)

plt.tight_layout()
plt.show()
