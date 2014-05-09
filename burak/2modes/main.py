"""
	Main file to produce results that are going to be included in the 2modes paper
	Not every calculation needs to be done at each run, data can be produced
	and stored locally, see the booleans in line 18
"""

import numpy as np
from scipy import interpolate, integrate
from scipy.misc import derivative
#Initiate plotting environment:
import matplotlib as mpl
from pylab import plot, xlabel, ylabel, show
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#Import twomode system module
import twomode 

#Booleans:
computeSolution = False
computePsect = False
computeArcLengths = False
plotPsect = False
plotRetmap = True

#Only relative equilibrium:
reqv = np.array([0.43996557973671596,
				 0,
				 -0.38626705952942764,
				 0.07020440068369532], float) 
#Compute reduced stability matrix for the relative equilibrium:
Ared = twomode.StabilityMatrixRed(reqv)
#Compute eigenvalues and eigenvectors:
w, v = np.linalg.eig(Ared)
#Pick the real part of unstable eigenvector as the Poincare section direction:
unstabledir = np.real(v[:,0])
#unstabledir = np.real(w[0]*v[:,0]+w[1]*v[:,1])
unstabledir = unstabledir/np.linalg.norm(unstabledir) #Normalize
#Auxiliary vector to define Poincare basis:
vaux3D = np.array([0,0,1], float)
unstabledir3D = twomode.four2three(unstabledir)
nhat3D = np.cross(vaux3D, unstabledir3D)
nhat = twomode.three2four(nhat3D)
vaux = twomode.three2four(vaux3D)
#Poincare Basis:
psbasis = np.array([unstabledir, vaux, nhat]).transpose()
			
if computeSolution:
	
	tf = 4000;
	dt = 0.001;
	epsilon = 1e-3;
	x0 = reqv+epsilon*unstabledir
	xphi0 = np.append(x0, 0)
	print xphi0
	
	t = np.linspace(0, tf, np.floor(tf/dt)+1)
	xphisol = twomode.intslice(xphi0, t)

	txphisol = np.concatenate((t.reshape(-1,1),xphisol), axis=1)
	#Create a figure window
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	#Modify axis colors:
	ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
	ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
	ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
			
	ax.plot(xphisol[:,0], 
	xphisol[:,2], 
	xphisol[:,3], linewidth=0.1, color='#3c5f96')
	
	np.savetxt('data/txphisol.dat', txphisol)
	
	plt.show()
	
if computePsect:
	import poincare
	if not('txphisol' in locals()):
		txphisol = np.loadtxt('data/txphisol.dat')
	ps = poincare.computeps(txphisol, reqv, nhat, 1)
	np.savetxt('data/ps.dat', ps)

if not('ps' in locals()):
	ps = np.loadtxt('data/ps.dat')

#Project onto poincare basis:
ps2D = np.dot(psbasis.transpose(), (ps[:,1:5]-reqv).transpose()).transpose()
ps2D = ps2D[:,0:2]
#Sort ps2D for interpolation:
sortingindices = np.argsort(ps2D[:,0])
ps2Dsorted = np.array(ps2D[sortingindices, :],float)

#Interpolate to psect:
print 'Interpolating the Poincare section'
tckps = interpolate.splrep(ps2Dsorted[:,0],ps2Dsorted[:,1])
dxintps = (ps2Dsorted[np.size(ps2Dsorted,0)-1,0] - ps2Dsorted[0,0])/100
xintps = np.arange(ps2Dsorted[0,0], 
					  ps2Dsorted[np.size(ps2Dsorted,0)-1,0]+dxintps, 
					  dxintps)
yintps = interpolate.splev(xintps, tckps)

def psarclength(x):
	"""
	Function to compute arclenghts in the Poincare section
	Ref: http://www.mathwords.com/a/arc_length_of_a_curve.htm 
	"""
	def ds(x):
		dypsect = derivative(interpolate.splev, x, dx=1e-6, args=(tckps,), order=5)
		return np.sqrt(1 + dypsect**2)
	quad, err = integrate.quad(ds, ps2Dsorted[0,0], x)	
	return quad

if computeArcLengths:
	print 'Computing arclengths corresponding to data'
	sn = np.array([psarclength(ps2D[i,0]) for i in range(np.size(ps2D,0))], float)
	snplus1 = sn[1:]
	sn = sn[0:-1]
	RetMapData = np.array([sn, snplus1], float).transpose()
	np.savetxt('data/RetMapData.dat', RetMapData)
	
if not('sn' in locals()):
	RetMapData = np.loadtxt('data/RetMapData.dat')
	sn = RetMapData[:,0]
	snplus1 = RetMapData[:,1]

print 'Interpolating to the return map'
isortRetMap = np.argsort(sn)
tckRetMap = interpolate.splrep(sn[isortRetMap],snplus1[isortRetMap], k=3)
dxintRetMap = (np.max(sn)-np.min(sn))/50000
xintRetMap = np.arange(np.min(sn), np.max(sn), dxintRetMap)
yintRetMap = interpolate.splev(xintRetMap, tckRetMap)

#Interpolation in two parts (Since the typica return map is discontinuous 
#in the first derivative):
print 'Interpolating to the return map in 2 parts'
snsorted = sn[isortRetMap]
snplus1sorted = snplus1[isortRetMap]
imax = np.argmax(snplus1sorted)
tck1 = interpolate.splrep(snsorted[0:imax+1], snplus1sorted[0:imax+1])
xint1 = np.linspace(snsorted[0], snsorted[imax], 50000)
yint1 = interpolate.splev(xint1, tck1)

tck2 = interpolate.splrep(snsorted[imax:len(snsorted)], snplus1sorted[imax:len(snsorted)])
xint2 = np.linspace(snsorted[imax], snsorted[len(snsorted)-1], 50000)
yint2 = interpolate.splev(xint2, tck2)

snint = np.append(xint1, xint2)
snplus1int = np.append(yint1, yint2)

#Return map function:
def retmap(sn):
	if sn < snsorted[imax]:
		snp1 = interpolate.splev(sn, tck1)
	else:
		snp1 = interpolate.splev(sn, tck2)
	return snp1
#nth return map function:
def retmapm(n, sn):
	snpn = retmap(sn)
	for i in range(n - 1):
		snpn = retmap(snpn)
	return	snpn


if plotPsect:
	plot(ps2D[:,0], ps2D[:,1], '.')
	plt.hold(True)
	plot(xintps, yintps)
	show()

if plotRetmap:
	plot(sn, snplus1, '.')
	plt.hold(True)
	plot(xintRetMap, yintRetMap)
	#plot(snint, snplus1int)
	show()
