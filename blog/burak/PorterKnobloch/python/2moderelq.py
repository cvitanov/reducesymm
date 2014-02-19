#
# invpolsolver.py
#
"""
Use odeint to solve differential equations defined by vinvpol in twomode.py
"""

from scipy.integrate import odeint
import twomode
import numpy as np
import sspsolver

import matplotlib as mpl
from pylab import figure, plot, xlabel, ylabel, grid, hold, legend, title, savefig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from subprocess import call

#Load parameters:
p = np.loadtxt('data/parameters.dat')

x01 = [0.43997, 0, -0.38627, 0.07020] #reqv

t1 = np.linspace(0,200,10000)

xsol1 = sspsolver.integrate(x01, p, t1)

np.savetxt('data/2moderpoBrown.dat', xsol1)

fig = plt.figure(figsize=(6,4.5))
ax = fig.gca(projection='3d')
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))

ax.plot(xsol1[:,0], xsol1[:,1], xsol1[:,2], linewidth=0.4, color='#a81800')

ax.set_xlabel('$x_1$', fontsize=14)
ax.set_ylabel('$y_2$', fontsize=14)
ax.set_zlabel('$x_2$', fontsize=14)

Nticks = 3

xticks = np.linspace(-0.3, 0.3, Nticks)
ax.set_xticks(xticks)
ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=8); # use LaTeX formatted labels

yticks = np.linspace(-0.3, 0.3, Nticks)
ax.set_yticks(yticks)
ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=8); # use LaTeX formatted labels

zticks = np.linspace(-0.25, 0.25, Nticks)
ax.set_zticks(zticks)
ax.set_zticklabels(["$%.1f$" % ztik for ztik in zticks], fontsize=8); # use LaTeX formatted labels

savefig('dasbuchfigs/2modereqv.png', bbox_inches='tight', dpi=75)
call(['convert', '-trim', "dasbuchfigs/2modereqv.png", "dasbuchfigs/2modereqv.png"])
call(["inkscape", "dasbuchfigs/2modereqv.png", "-E", "dasbuchfigs/2modereqv.eps"])

#plt.show()
