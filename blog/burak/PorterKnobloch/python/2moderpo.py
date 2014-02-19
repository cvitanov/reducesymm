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
x02 = [0.001, 0, 0.001, 0.001] #origin
x03 = [0.4525719,   0.0, 0.0509257, 0.0335428] #rpo
t1 = np.linspace(0,200,10000)
t2 = np.linspace(0,15,1500)
t3 = np.linspace(0, 3*3.6415120, 10000)

xsol1 = sspsolver.integrate(x01, p, t1)
xsol2 = sspsolver.integrate(x02, p, t2)
xsol3 = sspsolver.integrate(x03, p, t3)

x04 = xsol2[-1,:]
xsol4 = sspsolver.integrate(x04, p, t1/3)

np.savetxt('data/2moderpoBrown.dat', xsol1)
np.savetxt('data/2moderpoGreen.dat', xsol2)
np.savetxt('data/2moderpoMagenta.dat', xsol3)
np.savetxt('data/2moderpoBlue.dat', xsol4)

fig = plt.figure(figsize=(6,4.5))
ax = fig.gca(projection='3d')
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))

ax.plot(xsol1[:,0], xsol1[:,1], xsol1[:,2], linewidth=0.4, color='#a81800')
ax.hold(True)
ax.plot(xsol2[:,0], xsol2[:,1], xsol2[:,2], linewidth=1, color='#6a9c53')
ax.plot(xsol3[:,0], xsol3[:,1], xsol3[:,2], linewidth=0.8, color='#e500e5')
ax.plot(xsol4[:,0], xsol4[:,1], xsol4[:,2], linewidth=0.3, color='#3c5f96')
ax.set_xlabel('$x_1$', fontsize=14)
ax.set_ylabel('$y_2$', fontsize=14)
ax.set_zlabel('$x_2$', fontsize=14)

Nticks = 3

xticks = np.linspace(-1.5, 1.5, Nticks)
ax.set_xticks(xticks)
ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=8); # use LaTeX formatted labels

yticks = np.linspace(-1.5, 1.5, Nticks)
ax.set_yticks(yticks)
ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=8); # use LaTeX formatted labels

zticks = np.linspace(-1.5, 1.5, Nticks)
ax.set_zticks(zticks)
ax.set_zticklabels(["$%.1f$" % ztik for ztik in zticks], fontsize=8); # use LaTeX formatted labels

savefig('dasbuchfigs/py2moderpo.png', bbox_inches='tight', dpi=75)
call(['convert', '-trim', "dasbuchfigs/py2moderpo.png", "dasbuchfigs/py2moderpo.png"])
call(["inkscape", "dasbuchfigs/py2moderpo.png", "-E", "dasbuchfigs/py2moderpo.eps"])

#plt.show()
