from pylab import figure, grid, hold, legend, savefig
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 

import numpy as np

taux = np.loadtxt("data/solutionscaledtime.dat")
xsol = taux[:,1:5]

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

fig = plt.figure() #Create a figure instance
ax = fig.gca(projection = "3d") #Create an axis instance

ax.plot(xsol[:,0], xsol[:,2], xsol[:,3], lw=0.3)
ax.set_xlabel('$\hat{x}_1$', fontsize=18)
ax.set_ylabel('$\hat{x}_2$', fontsize=18)
ax.set_zlabel('$\hat{y}_2$', fontsize=18)
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
ax.view_init(25,35)
savefig("image/scaledtimesol.png", bbox_inches='tight', dpi=100)
