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
from scipy import interpolate, integrate
from scipy.misc import derivative
from scipy.optimize import newton

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
#sectp = np.array([0, 0, 0, 0], float)
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
	
	#arclength[i] = s

	#arclength[i] = farclength(psectsorted[i,0])

sn = np.zeros(np.size(psectsorted,0))
sn[sortedindices] = arclength
snplus1 = sn[1:len(sn)]
sn=sn[0:len(sn)-1]

snsortedindices = np.argsort(sn)
snsorted = sn[snsortedindices]
snplus1sorted = snplus1[snsortedindices]

#Interpolation in two parts:
imax = np.argmax(snplus1sorted)
tck1 = interpolate.splrep(snsorted[0:imax], snplus1sorted[0:imax])
xint1 = np.arange(snsorted[0], snsorted[imax], snsorted[imax]/100)
yint1 = interpolate.splev(xint1, tck1)

tck2 = interpolate.splrep(snsorted[imax:len(snsorted)-1], snplus1sorted[imax:len(snsorted)-1])
xint2 = np.arange(snsorted[imax], 
				  snsorted[len(snsorted)-1]+snsorted[imax]/100, 
				  snsorted[imax]/100)
yint2 = interpolate.splev(xint2, tck2)

#Return map function
def retmap(sn):
	
	if sn < snsorted[imax]:
		
		snp1 = interpolate.splev(sn, tck1)
		
	else:
		
		snp1 = interpolate.splev(sn, tck2)

	return snp1

#Higher order return map functions:
def retmap2(sn):
	
	return retmap(retmap(sn))

def retmap3(sn):
	
	return retmap(retmap(retmap(sn)))

def retmap4(sn):
	
	return retmap(retmap(retmap(retmap(sn))))
	
def retmap5(sn):
	
	return retmap(retmap(retmap(retmap(retmap(sn)))))


#Find RPOs:

#Function to solve:
def srpo1(sn):
	
	return retmap(sn) - sn

#Find the arclength corresponding to the rpo
srpo = newton(srpo1, 0.1)

#Guessing the initial point:
i0 = np.argmin(np.absolute(arclength - srpo))

if arclength[i0] < srpo:
	
	i1 = i0
	
else:
	
	i1 = i0-1

dspsect = arclength[i1+1] - arclength[i1]
dsrpo = srpo - arclength[i1]
dpsect = psectsorted[i1+1, 0:3] - psectsorted[i1, 0:3]

psrporel3D = psectsorted[i1, 0:3] + dpsect*(dsrpo/dspsect)


#One must solve the inverse arclength function to find the position of rpo
#on the psect:

#def arclength2psectrel(px, s):
	#'''
		#Function to be solved to get relative x coordinate of the rpo on
		#the Poincare section
	#'''
	##sfun, err = farclength(px)
	#sfun = farclength(px)
	
	#return sfun-s

##Guessing the initial point:
#i0 = np.argmin(np.absolute(arclength - srpo))
#p0x = psectsorted[i0,0]

#pxrpo = newton(arclength2psectrel, p0x, args=(srpo,))

#pyrpo = interpolate.splev(pxrpo, tckpsect)

def psrelproj2xhat(psrelproj3D):
	'''
		Takes the relative position projected on Poincare section basis 
		on the Poincare section and returns the position on the reduced 
		state space
	'''	
	psrel3D = np.dot(psectbasis.transpose(), psrelproj3D)
		
	ps3D = psrel3D + sectp3D
	
	ps = np.array([ps3D[0], 0, ps3D[1], ps3D[2]], float)
	
	return ps	

rpoxhat = psrelproj2xhat(psrporel3D)

print "srpo="
print srpo
print "Relative periodic orbit passes through:"
print "x1: %5.16f" % rpoxhat[0]
print "y1: %5.16f" % rpoxhat[1]
print "x2: %5.16f" % rpoxhat[2]
print "y2: %5.16f" % rpoxhat[3]

#Plotting:

lw = 1.5

#Plot Poincare section:
figure(1, figsize=(8, 8))

xlabel('$v_u$', fontsize=36)
ylabel('$e_{y_2}$', fontsize=36)

plot(psectprojected[:,0],psectprojected[:,1], '.', ms=5)
#plt.grid()
plt.hold(True)

#plot(xintpsect, yintpsect, c='k', linewidth=lw)

savefig('image/psectonslice.png', bbox_inches='tight', dpi=100)

srange = np.arange(np.min(sn), np.max(sn), (np.max(sn)-np.min(sn))/100)

sp1 = np.zeros(np.size(srange, 0))
sp2 = np.zeros(np.size(srange, 0))
sp3 = np.zeros(np.size(srange, 0))
sp4 = np.zeros(np.size(srange, 0))
sp5 = np.zeros(np.size(srange, 0))

for i in range(np.size(srange,0)):
	
	sp1[i] = retmap(srange[i])
	sp2[i] = retmap2(srange[i])
	sp3[i] = retmap3(srange[i])
	sp4[i] = retmap4(srange[i])
	sp5[i] = retmap5(srange[i])

#Plot return maps:
figure(2, figsize=(8, 8))

xlabel('$s_n$', fontsize=36)
ylabel('$s_{n+1}$', fontsize=36)

plot(snsorted,snplus1sorted, '.', ms=6)
#plt.grid()
plt.hold(True)

plot(np.arange(np.min(sn), np.max(sn), np.max(sn)/100), np.arange(np.min(sn), np.max(sn), np.max(sn)/100))
plot(srange, sp1, c='k')
#plot(xint1, yint1, c='k', linewidth=lw)
#plot(xint2, yint2, c='k', linewidth=lw)

mpl.rcParams['grid.linewidth'] = 2

savefig('image/retmaponslice.png', bbox_inches='tight', dpi=100)

#figure(3, figsize=(8, 8))

#plot(srange, sp2)
#plt.hold(True)
#plot(srange, srange, c='k')

#figure(4, figsize=(8, 8))

#plot(srange, sp3)
#plt.hold(True)
#plot(srange, srange, c='k')

#figure(5, figsize=(8, 8))

#plot(srange, sp4)
#plt.hold(True)
#plot(srange, srange, c='k')

#figure(6, figsize=(8, 8))

#plot(srange, sp5)
#plt.hold(True)
#plot(srange, srange, c='k')


#f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
#ax1.plot(np.arange(np.min(sn), np.max(sn), np.max(sn)/100), retmap2(np.arange(np.min(sn), np.max(sn), np.max(sn)/100)))
#ax1.hold(True)
#ax1.plot(np.arange(np.min(sn), np.max(sn), np.max(sn)/100), 
		 #np.arange(np.min(sn), np.max(sn), np.max(sn)/100))
#ax2.
#ax3.
#ax4.plot(x, 2 * y ** 2 - 1, color='r')

plt.tight_layout()
plt.show()
