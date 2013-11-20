#
#psectandretmap.py
#
"""
	Calculates the poincare section intersections using the previously computed
	solution on the slice.
	Plots the poincare section.
	Computes and plots the return map.
	This will produce correct results only if the xhatp=[1,0,0,0] slice 
	is used.
"""

import numpy as np
from scipy import interpolate

from pylab import figure, plot, xlabel, ylabel, grid, hold, legend, title, savefig
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

import twomode
import onslicesolver
import psectslice

pars = np.loadtxt('data/parameters.dat')

sectp = np.array([0.4399655797367152, 0, -0.38626706847930564, 0.0702043939917171], float)
sectp3D = np.array([sectp[0], sectp[2], sectp[3]],float)
vaux = np.array([0,0,1], float) #Auxiliary vector to decide the Poincare section direction
unstabledir = np.array([0.02884567,  0.99957642, -0.00386173], float) #Unstable direction
nhat3D = np.cross(vaux, unstabledir)
	
nhat = np.array([nhat3D[0],0,nhat3D[1],nhat3D[2]], float)
direction = 1 # set 1 or -1 to choose between directions of piercing the Poincare section

compute = 0

if compute:	
	xhat = np.loadtxt('data/solutiononslice.dat') #Solution on slice
	#compute Poincare section:
	ps=psectslice.computeps(xhat, sectp, nhat, direction,  p=pars)
else:
	ps = np.loadtxt('data/psectslice.dat')

#Write non-zero values on a (:,3) array:
ps3D = np.array([ps[:,1], ps[:,3], ps[:,4]],float).transpose()
#Calculate positions relative to the section point:
psrel = ps3D -  sectp3D
#Generate the matrix with the basis of Poincare section hyperplane:
psectbasis = np.array([unstabledir,vaux,nhat3D],float)
#Project Poincare section onto section basis:
psectprojected = np.dot(psectbasis,psrel.transpose()).transpose()

sortedindices = np.argsort(psectprojected[:,0])
psectsorted = np.array(psectprojected[sortedindices, :],float)

arclength = np.zeros(np.size(psectsorted,0))

for i in range(1,np.size(psectsorted,0)):
	arclength[i] = arclength[i-1] + np.linalg.norm(psectsorted[i,:]-psectsorted[i-1,:])

sn = np.zeros(np.size(psectsorted,0))
sn[sortedindices] = arclength
snplus1 = sn[1:len(sn)]
sn=sn[0:len(sn)-1]

snsortedindices = np.argsort(sn)
snsorted = sn[snsortedindices]
snplus1sorted = snplus1[snsortedindices]

#Interpolation in two parts:
imax = np.argmax(snplus1sorted)
tck1 = interpolate.splrep(snsorted[0:imax+1], snplus1sorted[0:imax+1])
xint1 = np.arange(snsorted[0], snsorted[imax+1], snsorted[imax]/100)
yint1 = interpolate.splev(xint1, tck1)

tck2 = interpolate.splrep(snsorted[imax:len(snsorted)-1], snplus1sorted[imax:len(snsorted)-1])
xint2 = np.arange(snsorted[imax], snsorted[len(snsorted)-1], snsorted[imax]/100)
yint2 = interpolate.splev(xint2, tck2)

#Plot Poincare section:
figure(1, figsize=(8, 8))

xlabel('$v_u$', fontsize=36)
ylabel('$e_{y_2}$', fontsize=36)

plot(psectprojected[:,0],psectprojected[:,1], '.', ms=5)
plt.grid()

savefig('image/psectonslice.png', bbox_inches='tight', dpi=100)

#Plot return map:
figure(2, figsize=(8, 8))

xlabel('$s_n$', fontsize=36)
ylabel('$s_{n+1}$', fontsize=36)

plot(snsorted,snplus1sorted, '.', ms=6)
plt.grid()
plt.hold(True)

lw = 1.5
plot(np.arange(np.min(sn), np.max(sn), np.max(sn)/100), np.arange(np.min(sn), np.max(sn), np.max(sn)/100))
plot(xint1, yint1, c='k', linewidth=lw)
plot(xint2, yint2, c='k', linewidth=lw)

mpl.rcParams['grid.linewidth'] = 2

savefig('image/retmaponslice.png', bbox_inches='tight', dpi=100)

plt.tight_layout()
plt.show()
