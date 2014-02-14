#
# plotfull.py
#

from __future__ import unicode_literals
import KS
import numpy as np
from scipy.linalg import expm

tx = np.loadtxt("data/solutionscaledtime.dat") #Load the scaled time solution

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

ax.plot(tx[:,1], tx[:,3], tx[:,4], lw=0.5)
#ax.plot(tx[:,10], tx[:,11], tx[:,12], lw=0.5)
ax.set_xlabel('$\hat{x}_1$', fontsize=18)
ax.set_ylabel('$\hat{x}_2$', fontsize=18)
ax.set_zlabel('$\hat{y}_2$', fontsize=18)
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
ax.view_init(20,-65)
savefig("image/sspRedScaledTime.png", bbox_inches='tight', dpi=100)

plt.tight_layout()
plt.show()
