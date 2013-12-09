#
# rpoplotter.py
#
"""
Plot relative periodic orbits
"""
#Library imports:
from __future__ import unicode_literals
from numpy import loadtxt
from pylab import figure, plot, xlabel, ylabel, grid, hold, legend, title, savefig
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
#User module imports:
import onslicesolver

#Load parameters:
pars = np.loadtxt('data/parameters.dat')

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

#Read integration times:
Trpo = np.loadtxt('data/RPOT.dat')

#Read RPOs:
xrpo = np.loadtxt('data/50rpos.dat')
#xrpo = np.loadtxt('data/rposxhat.dat')

#How many RPO's to plot:
nplot = 2

abserr = 1.0e-14
relerr = 1.0e-13

#Start plotting:
fig = plt.figure()
ax = fig.gca(projection='3d')

lw = 1
axfsize = 36

ax.set_xlabel('$x_1$', fontsize=axfsize)
ax.set_ylabel('$y_1$', fontsize=axfsize)
ax.set_zlabel('$y_2$', fontsize=axfsize)	

ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
ax.view_init(30,30)

ax.hold(True)
	
for i in range(nplot):
	
	stoptime = Trpo[i]
	numpoints = int(stoptime/0.001 + 1.0)
	
	t = np.linspace(0,stoptime,numpoints)
	
	xphi0=[xrpo[i, 0], 0, xrpo[i, 1], xrpo[i, 2], 0]
	#xphi0=[xrpo[i, 0], xrpo[i, 1], xrpo[i, 2], xrpo[i, 3], 0]
	
	xsol = onslicesolver.integrate(xphi0, pars, t, abserror=abserr, relerror=relerr)
	
	if i%nplot == 0:
		color = 'b'
	elif i%nplot == 1:
		color = 'r'
	elif i%nplot == 2:
		color = 'g'
	elif i%nplot == 3:
		color = 'k'
		
	ax.plot(xsol[:,0], xsol[:,2], xsol[:,3], linewidth=lw , c=color)
	
	print xsol[0,:] - xsol[numpoints-1,:]

savefig('image/rpos.png', bbox_inches='tight', dpi=150)
	
plt.tight_layout()
plt.show()
