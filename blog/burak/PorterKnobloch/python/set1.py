#
# ssp3dplotter.py
#
"""
Plot the solution generated by invpolsolver.py
"""

from __future__ import unicode_literals
from numpy import loadtxt
from pylab import figure, plot, xlabel, ylabel, grid, hold, legend, title, savefig
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

t, x1, y1, x2, y2 = loadtxt('data/sspsolution.dat', unpack=True)
t, u, v, w, q = loadtxt('data/invpolsolution.dat', unpack=True)
t, x1hat, y1hat, x2hat, y2hat, phi = loadtxt('data/solutiononslice.dat', unpack=True)

#Define linewidth:
lw = 0.6
axfsize = 36

fig = plt.figure()

ax = fig.gca(projection='3d')

ax.plot(x1[10000:50000], y1[10000:50000], y2[10000:50000], linewidth=lw)
ax.set_xlabel('$x_1$', fontsize=axfsize)
ax.set_ylabel('$y_1$', fontsize=axfsize)
ax.set_zlabel('$y_2$', fontsize=axfsize)
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
ax.view_init(30,30)

ax.hold(True)
ax.plot(x1[0:10000], y1[0:10000], y2[0:10000], linewidth=lw, c='r')
savefig('image/Set1ssp1.png', bbox_inches='tight', dpi=100)

fig2 = plt.figure()

ax2 = fig2.gca(projection='3d')

ax2.plot(x1[10000:50000], x2[10000:50000], y2[10000:50000], linewidth=lw)
ax2.set_xlabel('$x_1$', fontsize=axfsize)
ax2.set_ylabel('$x_2$', fontsize=axfsize)
ax2.set_zlabel('$y_2$', fontsize=axfsize)
ax2.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax2.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax2.w_zaxis.set_pane_color((1, 1, 1, 1.0))
ax2.view_init(30,30)

ax2.hold(True)
ax2.plot(x1[0:10000], x2[0:10000], y2[0:10000], linewidth=lw, c='r')
savefig('image/Set1ssp2.png', bbox_inches='tight', dpi=100)

fig3 = plt.figure()
ax3 = fig3.gca(projection='3d')

ax3.plot(u[10000:50000], v[10000:50000], w[10000:50000], linewidth=0.8*lw)
#ax.grid(False)
ax3.set_xlabel('$u$', fontsize=axfsize)
ax3.set_ylabel('$v$', fontsize=axfsize)
ax3.set_zlabel('$w$', fontsize=axfsize)
ax3.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax3.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax3.w_zaxis.set_pane_color((1, 1, 1, 1.0))
ax3.view_init(15,30)

ax3.hold(True)
ax3.plot(u[0:10000], v[0:10000], w[0:10000], linewidth=5*lw, c='r')
savefig('image/Set1invpol.png', bbox_inches='tight', dpi=100)

fig4 = plt.figure()
ax4 = fig4.gca(projection='3d')

ax4.plot(x1hat[10000:50000], x2hat[10000:50000], y2hat[10000:50000], linewidth=0.5*lw)
ax4.set_xlabel('$\hat{x}_1$', fontsize=axfsize)
ax4.set_ylabel('$\hat{x}_2$', fontsize=axfsize)
ax4.set_zlabel('$\hat{y}_2$', fontsize=axfsize)
ax4.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax4.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax4.w_zaxis.set_pane_color((1, 1, 1, 1.0))
ax4.view_init(20, 30)

ax4.hold(True)
ax4.plot(x1hat[0:10000], x2hat[0:10000], y2hat[0:10000], linewidth=5*lw, c='r')
savefig('image/Set1sspred.png', bbox_inches='tight', dpi=100)
	
#plt.tight_layout()
#plt.show()
