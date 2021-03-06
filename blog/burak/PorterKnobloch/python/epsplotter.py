#
# ssp3dplotter.py
#
"""
Plot the solution generated by invpolsolver.py
"""

from __future__ import unicode_literals
from numpy import loadtxt
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
#mpl.use('ps')
from pylab import figure, plot, xlabel, ylabel, grid, hold, legend, title, savefig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

from subprocess import call

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

t, x1, y1, x2, y2 = loadtxt('data/sspsolution.dat', unpack=True)
t, x1hat, y1hat, x2hat, y2hat, phi = loadtxt('data/solutiononslice.dat', unpack=True)


fig = plt.figure(figsize=(4,3), dpi=10)
ax = fig.gca(projection='3d')

points = range(int(len(x1))-30000, int(len(x1)))

ax.plot(x1[points], x2[points], y2[points], lw=0.0001, markevery = 5)
ax.set_xlabel('$x_1$', fontsize=10)
ax.set_ylabel('$x_2$', fontsize=10)
ax.set_zlabel('$y_2$', fontsize=10)
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
ax.view_init(30,30)

Nticks = 5

xticks = np.linspace(np.min(x1), np.max(x1), Nticks)
ax.set_xticks(xticks)
ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=8); # use LaTeX formatted labels

yticks = np.linspace(np.min(x2), np.max(x2), Nticks)
ax.set_yticks(yticks)
ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=8); # use LaTeX formatted labels

zticks = np.linspace(np.min(y2), np.max(y2), Nticks)
ax.set_zticks(zticks)
ax.set_zticklabels(["$%.1f$" % ztik for ztik in zticks], fontsize=8); # use LaTeX formatted labels

savefig('image/ssp.eps', bbox_inches='tight')

call(["eps2eps", "image/ssp.eps", "image/ssp2.eps"])

fig.clf()

points = range(int(len(x1))-90000, int(len(x1)-60000))

ax = fig.gca(projection='3d')
ax.plot(x1hat[points], x2hat[points], y2hat[points], linewidth=0.3)
ax.set_xlabel('$x_1$', fontsize=10)
ax.set_ylabel('$x_2$', fontsize=10)
ax.set_zlabel('$y_2$', fontsize=10)
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
ax.view_init(30,30)

Nticks = 5

xticks = np.linspace(np.min(x1hat), np.max(x1hat), Nticks)
ax.set_xticks(xticks)
ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=8); # use LaTeX formatted labels

yticks = np.linspace(np.min(x2hat), np.max(x2hat), Nticks)
ax.set_yticks(yticks)
ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=8); # use LaTeX formatted labels

zticks = np.linspace(np.min(y2hat), np.max(y2hat), Nticks)
ax.set_zticks(zticks)
ax.set_zticklabels(["$%.1f$" % ztik for ztik in zticks], fontsize=8); # use LaTeX formatted labels

savefig('image/sspred.eps', bbox_inches='tight')

call(["eps2eps", "image/sspred.eps", "image/sspred2.eps"])

	
#plt.tight_layout()
#plt.show()
