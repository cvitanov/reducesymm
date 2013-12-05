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
import rpo

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

#Interpolate to psect:
tckpsect = interpolate.splrep(psectsorted[:,0],psectsorted[:,1])
dxintpsect = (psectsorted[np.size(psectsorted,0)-1,0] - psectsorted[0,0])/100
xintpsect = np.arange(psectsorted[0,0], 
					  psectsorted[np.size(psectsorted,0)-1,0]+dxintpsect, 
					  dxintpsect)
yintpsect = interpolate.splev(xintpsect, tckpsect)

#Functions to compute arclengths:
#A ref: http://www.mathwords.com/a/arc_length_of_a_curve.htm
def dypsect(x):
	
	return derivative(interpolate.splev, x, dx=1e-6, args=(tckpsect,), order=5)

def ds(x):
	
	return np.sqrt(1 + dypsect(x)**2)

def farclength(x):
	
	simpson = 0
	
	if simpson:
	
		dx = (psectsorted[1,0]-psectsorted[0,0])/10
		xrang = np.arange(psectsorted[0,0]-2*dx, x+2*dx, dx)
		yint = interpolate.splev(xrang, tckpsect)
		
		#http://en.wikipedia.org/wiki/Five-point_stencil#First_derivative
		ny = np.size(yint,0)
		dydx = (-yint[2+2:ny] + 8*yint[2+1:ny-1] - 8*yint[2-1:ny-3] + yint[0:ny-4])/(12*dx)
		dsdx = np.sqrt(1+dydx**2)
		
		simps = integrate.simps(dsdx, xrang[2:ny-2])
	
		return simps
	
	else:
	
		quad, err = integrate.quad(ds, psectsorted[0,0], x)	
	
		return quad 		

arclength = np.zeros(np.size(psectsorted,0))

for i in range(1,np.size(psectsorted,0)):
	
	#arclength[i] = arclength[i-1] + np.linalg.norm(psectsorted[i,:]-psectsorted[i-1,:])
	#s, err = farclength(psectsorted[i,0])
	#arclength[i] = s
	
	arclength[i] = farclength(psectsorted[i,0])

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
xint1 = np.arange(snsorted[0], snsorted[imax], snsorted[imax]/100)
yint1 = interpolate.splev(xint1, tck1)

tck2 = interpolate.splrep(snsorted[imax:len(snsorted)], snplus1sorted[imax:len(snsorted)])
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

def retmapn(n, sn):
	snpn = retmap(sn)
	for i in range(n - 1):
		snpn = retmap(snpn)
	return	snpn
	
#Find RPOs:

#Function to solve:
def srpo1(sn):
	return retmap(sn) - sn
def srpon(n,sn):
	return retmapn(n,sn) - sn

#One must solve the inverse arclength function to find the position of rpo
#on the psect:

def arclength2psectrel(px, s):
	'''
		Function to be solved to get relative x coordinate of the rpo on
		the Poincare section
	'''
	#sfun, err = farclength(px)
	sfun = farclength(px)
	
	return sfun-s


def psrelproj2xhat(pxrel, pyrel):
	'''
		Takes the relative position projected on Poincare section basis 
		on the Poincare section and returns the position on the reduced 
		state space
	'''
	psrelproj3D = np.array([pxrel, pyrel, 0], float)
	
	psrel3D = np.dot(psectbasis.transpose(), psrelproj3D)
		
	ps3D = psrel3D + sectp3D
	
	ps = np.array([ps3D[0], 0, ps3D[1], ps3D[2]], float)
	
	return ps	


#Number of RPOs to look for:
nrpo = 10;
nretmap = 1;

#Found PROs:
frpo = 0;

#Dummy assignment:
scandidates = [[0, 0]]

srange = np.arange(np.min(sn), np.max(sn), (np.max(sn)-np.min(sn))/1000)

#Find arclengths of RPO candidates:

while frpo < nrpo:
	sp = np.array([retmapn(nretmap, sn) for sn in srange])
	#Look for rpo candidates in the retmap:
	retmapminid = sp-srange
	for i in range(len(retmapminid) - 1):
		if retmapminid[i]*retmapminid[i+1] < 0:
			#If there is a zero crossing, take the arclength as a candidate for
			#being a fixed point of the returnmap:
			scc = srange[i]
			#Is this a new candidate?
			new = 1
			for j in range(len(scandidates)):
				if scc == scandidates[j][0]:
					new = 0
			if new:
				scandidates = np.append(scandidates, [[scc, nretmap]], axis=0)
				frpo = frpo + 1
	nretmap = nretmap + 1

scandidates = scandidates[1:len(scandidates), :]

#Find psect intersections corresponding to the arclengths of RPO candidates:

rpocandidatesxt = np.zeros([nrpo, 5])
savetxt('data/rpocandidates.dat', rpocandidatesxt)

for i in range(nrpo):
	#Guessing the initial point:
	i0 = np.argmin(np.absolute(arclength - scandidates[i, 0]))
	p0x = psectsorted[i0,0]
	
	pxrpo = newton(arclength2psectrel, p0x, args=(scandidates[i, 0],))
	pyrpo = interpolate.splev(pxrpo, tckpsect)
	rpoxhat = psrelproj2xhat(pxrpo, pyrpo)
	
	inotsorted = sortedindices[i0]
	trpo = ps[inotsorted + scandidates[i,1] ,0] - ps[inotsorted ,0]
	
	rpoxhatt = np.append(rpoxhat, trpo)
	
	rpocandidatesxt[i,:] = rpoxhatt 

rposxt = np.zeros([nrpo, 5]) 

for i in range(nrpo):
	rposxt[i,:] = rpo.findrpo(rpocandidatesxt[i,:])

np.savetxt('data/rpo.dat', rposxt)

#Plotting:

lw = 1.5

#Plot Poincare section:
figure(1, figsize=(8, 8))

xlabel('$v_u$', fontsize=36)
ylabel('$e_{y_2}$', fontsize=36)

plot(psectprojected[:,0],psectprojected[:,1], '.', ms=5)
#plt.grid()
plt.hold(True)

plot(xintpsect, yintpsect, c='k', linewidth=lw)

savefig('image/psectonslice.png', bbox_inches='tight', dpi=100)

#Plot return maps:
figure(2, figsize=(8, 8))

xlabel('$s_n$', fontsize=36)
ylabel('$s_{n+1}$', fontsize=36)

plot(snsorted,snplus1sorted, '.', ms=6)
#plt.grid()

plt.hold(True)

sp1 = np.array([retmapn(1, sn) for sn in srange])

plot(srange, srange)
plot(srange, sp1, c='k')
#plot(xint1, yint1, c='k', linewidth=lw)
#plot(xint2, yint2, c='k', linewidth=lw)

#figure(3, figsize=(8, 8))
#xlabel('$s_n$', fontsize=36)
#ylabel('$s_{n+2}$', fontsize=36)

#plt.hold(True)

#plot(np.arange(np.min(sn), np.max(sn), np.max(sn)/100), np.arange(np.min(sn), np.max(sn), np.max(sn)/100))
#plot(srange, sp1, c='k')

#mpl.rcParams['grid.linewidth'] = 2

savefig('image/retmaponslice.png', bbox_inches='tight', dpi=100)

plt.tight_layout()
plt.show()
