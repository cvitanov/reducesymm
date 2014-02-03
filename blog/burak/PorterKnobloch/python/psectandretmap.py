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
	
	This program also compute relative periodic orbits and symbolic dynamics.
	It got rather messy, but you can choose what to do or not to do by assigning
	binaries below.
"""

import numpy as np
from scipy import interpolate, integrate
from scipy.misc import derivative
from scipy.optimize import newton, fsolve, root

from pylab import figure, plot, xlabel, ylabel, grid, hold, legend, title, savefig
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

import twomode
import onslicesolver
import sspsolver
import movingframes
import psectslice
import rpo

#What to do what not to do:
computeps = False
computerpo = True
computesymbdyn = True
plotpsandretmap = True
solveFPO = False

m = 6

pars = np.loadtxt('data/parameters.dat')

sectp = np.array([0.4399655797367152, 0, -0.38626706847930564, 0.0702043939917171], float)
sectp3D = np.array([sectp[0], sectp[2], sectp[3]],float)
vaux = np.array([0,0,1], float) #Auxiliary vector to decide the Poincare section direction
unstabledir = np.array([0.02884567,  0.99957642, -0.00386173], float) #Unstable direction
nhat3D = np.cross(vaux, unstabledir)
	
nhat = np.array([nhat3D[0],0,nhat3D[1],nhat3D[2]], float)
direction = 1 # set 1 or -1 to choose between directions of piercing the Poincare section

if computeps:
	print 'Computing Poincare section'	
	xhat = np.loadtxt('data/solutiononslice.dat') #Solution on slice
	#compute Poincare section:
	ps=psectslice.computeps(xhat, sectp, nhat, direction,  p=pars)
	#Save it:
	np.savetxt('data/psectslice.dat', ps)
else:
	print 'Loading Poincare section data'
	ps = np.loadtxt('data/psectslice.dat')

#Generate the matrix with the basis of Poincare section hyperplane:
psbasis = np.array([unstabledir,vaux,nhat3D],float)	

def psxhat2ps2D(ps, psectbasis = psbasis):
	"""
	Inputs Poincare section coordinates within the slice and outputs 2D
	projection onto the Poincare section basis
	"""
	#Write non-zero values on a (:,3) array:
	if np.shape(ps)==(5,):
		ps3D = np.array([ps[1], ps[3], ps[4]],float)
	else:
		ps3D = np.array([ps[:,1], ps[:,3], ps[:,4]],float).transpose()
		
	#Calculate positions relative to the section point:
	psrel = ps3D -  sectp3D
	#Project Poincare section onto section basis:
	psectprojected = np.dot(psectbasis,psrel.transpose()).transpose()
	
	return psectprojected

psectprojected = psxhat2ps2D(ps, psbasis)

sortedindices = np.argsort(psectprojected[:,0])
psectsorted = np.array(psectprojected[sortedindices, :],float)

#Interpolate to psect:
print 'Interpolating the Poincare section'
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
	
	simpson = False
	
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

#Compute arclength values:
print 'Computing arclengths corresponding to data'
arclength = np.array([farclength(psectsorted[i,0]) for i in range(np.size(psectsorted,0))], float)

sn = np.zeros(np.size(psectsorted,0)) #Dummy assignment to sn
sn[sortedindices] = arclength
snplus1 = sn[1:len(sn)] 
sn=sn[0:len(sn)-1]

#Sort the return map data increasing in the value of sn, this is necessary
#for interpolating to the data:
snsortedindices = np.argsort(sn)
snsorted = sn[snsortedindices]
snplus1sorted = snplus1[snsortedindices]

#Interpolation in two parts (Since the typica return map is discontinuous 
#in the first derivative):
print 'Interpolating to the return map'
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
	
#Find RPOs:

#Look upto the mth return map to find the periodic orbit candidates:

scandidates = np.zeros([1,2]) #dummy matrix to hold po candidate arclengths
smin = np.min(snsorted)
smax = np.max(snsorted)

for i in range(m):
	#Define the function zeros of which would correspond to the periodic 
	#orbit arclengths:
	def fpo(s):
		po = retmapm(i+1, s) - s
		return po
	
	#print "retmap %i" %(i+1)
	
	fpoevo=0
	for s0 in np.arange(smin, smax, (smax-smin)/50000):
		fpoev = fpo(s0)
		
		if fpoev * fpoevo < 0: #If there is a zero-crossing, look for the root:	
			sc = newton(fpo, s0)
			#print "sc = %f" %sc
			
			newcandidate = 1
			for j in range(np.size(scandidates,0)):
				#Discard if found candidate is previously found:
				if np.abs(scandidates[j,1] - sc)<1e-9: 
					newcandidate=0
				
			if newcandidate:
				scandidates = np.append(scandidates, np.array([[i+1, sc]], float), axis=0)
			
		fpoevo = fpoev

scandidates = scandidates[1:np.size(scandidates,0), :]
print "Periodic orbit candidate arclengths:"
print scandidates

#Read the maximum value of the return map to determine the point of the 
#binary partition:

imaxint = np.argmax(snplus1int)
sboundary = snint[imaxint]

#Dummy itinerary array:
scitineraries = np.array([''])
#Determine itineraries:
for i in range(np.size(scandidates,0)):
	
	scandidate = scandidates[i,1]
	
	if scandidate < sboundary:
		itinerary = "0"
	else:
		itinerary = "1"
		
	m = int(scandidates[i,0])	
	
	for i in range(1,m):
		
		scandidate = retmap(scandidate) 
		if scandidate < sboundary:
			itinerary = itinerary+"0"
		else:
			itinerary = itinerary+"1"
	
	scitineraries = np.append(scitineraries, np.array([itinerary]))

scitineraries = scitineraries[1:]

print "itineraries:"
print scitineraries

#Dummy assignments:
group = np.zeros(len(scitineraries), int)
position = np.zeros(len(scitineraries), int)
#Group candidates into multiple shooting candidates according to their itineraries:
for i in range(np.size(scandidates,0)):
		
	"ith itinerary:"
	
	iti = scitineraries[i]
	#print "range(length of iti - 1)"
	#print range(len(iti)-1)

	#Indices of mth order periodic orbit candidates:
	imth = np.argwhere(scandidates[:,0]==scandidates[i,0])
	
	if len(imth)==1:
		group[i] = int(1)
		position[i] = int(1)
						
	if group[i] == 0:
		group[i] = int(np.max(group))+1 #ith candidate starts a new group
		position[i] = int(1)
		
		for shift in range(-len(iti)+1,0):
		
			for j in imth:
			#print j
			
				if scitineraries[j]  == (iti[shift:] + iti[0:len(iti)+shift]):
					group[j] = group[i] #jth candidate belong to the same group
					position[j] = shift + len(iti)+1
									
print "group:"		
print group		
print "position:"
print position

together = np.append(np.array([scitineraries]).transpose(), np.array([group]).transpose(), axis=1)
together = np.append(together, np.array([position]).transpose(), axis=1)
together = np.append(together, np.array([scandidates[:,1]], str).transpose(), axis=1)

print "Itinerary | Group | Position in the group"
print together

raw_input("Press Enter to continue...")

def fs2ps2D(px, s):
	'''
		Function to be solved to get relative x coordinate of the rpo on
		the Poincare section
	'''
	#sfun, err = farclength(px)
	sfun = farclength(px)
	
	return sfun-s

#From arclength to projected Poincare section coordinates:
def s2ps2D(s):
	
	#Guessing the initial point:
	i0 = np.argmin(np.absolute(arclength - s))
	p0x = psectsorted[i0,0]
	
	px = newton(fs2ps2D, p0x, args=(s,))
	py = interpolate.splev(px, tckpsect)
	
	ps2D = np.array([px, py])
	
	return ps2D

ps2Dcandidates = np.array([s2ps2D(s) for s in scandidates[:,1]] )
print "ps2Dcandidates:"
print ps2Dcandidates

def ps2D2psxhat(ps2D):
		'''
			Takes the relative position projected on Poincare section basis 
			on the Poincare section and returns the position on the reduced 
			state space
		'''
		psrelproj3D = np.array([ps2D[0], ps2D[1], 0], float)
		
		psrel3D = np.dot(psbasis.transpose(), psrelproj3D)
			
		ps3D = psrel3D + sectp3D
		
		ps = np.array([ps3D[0], 0, ps3D[1], ps3D[2]], float)
		
		return ps

pscandidates = np.array([ps2D2psxhat(ps2D) for ps2D in ps2Dcandidates] )
print "pscandidates"
print pscandidates

def timeofflight(xhat0):
	"""
	Computes time of flight for xhat0 on the Poincare section to come back onto
	the section for a second time 
	"""
	#Find the point on the computed poincare section, closest to the x0
	imin = np.argmin(np.linalg.norm(ps[:,1:5]-xhat0, axis=1))
	#Take its time of flight as
	if imin < np.size(ps,0)-1: 
		Tapproximate = ps[imin+1,0]-ps[imin,0]
	else:
		Tapproximate = ps[imin,0]-ps[imin-1,0]
	print "Tapproximate:"
	print Tapproximate
	#Integrate for a little bit longer than the approximated integration time:
	stoptime = 1.2*Tapproximate
	numpoints = int(stoptime/0.01)
	#Integration time array:
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	xhatphi0 = np.append(xhat0, np.array([0], float))
	xhatsol = onslicesolver.integrate(xhatphi0, pars, t, abserror=1.0e-14, relerror=1.0e-12)
	tx = np.append(np.array([t], float).transpose(), xhatsol, axis=1)
	#Compute Poincare section:
	psreturn=psectslice.computeps(tx, sectp, nhat, direction,  p=pars)
	print "psreturn:"
	print psreturn
	#Take the time nearest to the approximated time. This is due to the
	#fact that the array ps sometimes includes the initial point and sometimes
	#does not, hence we are not always sure the position of the first return.
	itof = np.argmin(np.abs(psreturn[:,0]-Tapproximate))
	tof = psreturn[itof,0]
	
	return tof

def phireturn(xhat0, tof):
	"""
	Computes phi for xhat0 on the Poincare section to come back onto
	the section for a second time 
	"""

	#stoptime = tof
	#numpoints = 2
	##Integration time array:
	#t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	#xhatphi0 = np.append(xhat0, np.array([0], float))
	#xhatsol = onslicesolver.integrate(xhatphi0, pars, t, abserror=1.0e-14, relerror=1.0e-12)
	#print "xhatsol"
	#print xhatsol
	#phi = xhatsol[1,4]
	

	stoptime = tof
	numpoints = 2
	#Integration time array:
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	
	xsol = sspsolver.integrate(xhat0, pars, t, abserror=1.0e-14, relerror=1.0e-12)
	phi = movingframes.x2xhatphi(xsol[1,:])[0,4]
	
	#print "xhatsol"
	#print xhatsol
	#phi = xhatsol[1,4]
	
	return phi

TOF = np.array([timeofflight(x0) for x0 in pscandidates], float)
print "time of flights:"
print TOF
phiret = np.array([phireturn(pscandidates[i,:], TOF[i]) for i in range(np.size(TOF,0))], float)
print "phiret:"
print phiret


np.savetxt('data/group.dat', group)
np.savetxt('data/itineraries.dat', scitineraries, fmt="%s")
np.savetxt('data/position.dat', position)
np.savetxt('data/tof0.dat',TOF)
np.savetxt('data/phi0.dat', phiret)
np.savetxt('data/xrpo0.dat', pscandidates)


def FPO(xn, tof):
	"""
	Vector input - Vector output function to solve for to get fixed points 
	or n-cycles of the Poincare return map. 
	xn = (n*(d-1))-dimensional array with n-cycle candidates at each row.
	xn = [x1_1
		  x1_2
		  x2_1
		  x2_2
		  x3_1
		  x3_2
		  ...
			]
	Where in xi_j first index denotes ith element of cycle and the second 
	index denotes jth coordinate of ith element.
	tof = n - dimensional vector to hold return time for each candidate
	"""
	#Reshape xn into a form such that at each row there are coordinates of
	#one element of the cycle in collumns. Namely,
	#xn = [[x1_1, x1_2], 
	#	   [x2_1, x2_2], 
	#	   [x3_1, x3_2], ...]
	xn = xn.reshape(np.size(tof), 2)
	#First compute P(x_i) and write as a matrix (n x d-1):
	Px = np.zeros(np.shape(xn)) #dummy assignment
	#print "Px"
	#print Px
	for i in range(np.size(xn,0)):
		#3D initial point to integrate
		xhat0=ps2D2psxhat(xn[i])
		
		pret = np.array([], float)
		factor = 1.2
		
		while np.size(pret)==0:
			#Integrate for a little bit longer than the time of flight:
			stoptime = factor*tof[i]
			#print "stoptime="
			#print stoptime
			numpoints = int(stoptime/0.01)
			#Integration time array:
			t = [stoptime * float(j) / (numpoints - 1) for j in range(numpoints)]
			xhatphi0 = np.append(xhat0, np.array([0], float))
			xhatsol = onslicesolver.integrate(xhatphi0, pars, t, abserror=1.0e-14, relerror=1.0e-12)
			tx = np.append(np.array([t], float).transpose(), xhatsol, axis=1)
			#Compute Poincare section:
			psreturn=psectslice.computeps(tx, sectp, nhat, direction,  p=pars)
			#Take the time nearest to the approximated time. This is due to the
			#fact that the array ps sometimes includes the initial point and sometimes
			#does not, hence we are not always sure the position of the first return.
			#print "xsol:"
			#print xsol
			#print "psreturn:"
			#print psreturn
			if np.size(psreturn,0) != 0:
				itof = np.argmin(np.abs(psreturn[:,0])-tof[i])
				pret = psreturn[itof, 0:]			
			else:
				factor = factor**2
				if factor > 50:
					print "something went wrong for:"
					print "xhat0 = "
					print xhat0
					print "Time of flight:"
					print tof[i]
					break 
					
		ps2D = psxhat2ps2D(pret)[0:2]
		
		#print "ps2D"
		#print ps2D
		
		Px[i,:] = ps2D
		
	Pxshifted = np.append(np.array([Px[np.size(Px,0)-1,:]]), Px[0:np.size(Px,0)-1,:], axis=0)	
	
	fun = xn - Pxshifted
	#Convert back to a vector function
	fun = fun.reshape(np.size(fun))
	
	return fun

if solveFPO:

	#Find periodic orbits from the cycles on the Poincare return map:
	x2Dpo = np.zeros(np.shape(ps2Dcandidates))
	convergence = np.array([''])
	for i in range(1,int(np.max(group)) +1):
		print "Group no:"
		print i
		#Get corresponding indices for the group i:
		indices = np.argwhere(group == i)
		#Reshape them in a 1D array:
		indices = indices.reshape(np.size(indices))
		#Get coordinates of group i
		xn = ps2Dcandidates[indices,:]
		#Get the positions in group i
		posn = position[indices]
		#Get time of flights for group i:
		tof = TOF[indices]
		#Sort tof and xn with respect to posn:
		xn = xn[posn-1,:]
		tof = tof[posn-1]
		#Turn xn to a 1xn(d-1) vector:	
		xn = xn.reshape(np.size(xn))
		#using fsolve:
		fsolveoutput = fsolve(FPO, xn, args=(tof), full_output=1, xtol=1e-10)
		print fsolveoutput
		xpo = fsolveoutput[0]
		convergence = np.append(convergence, np.array([fsolveoutput[3]]))
		print fsolveoutput[3] 
	
		#using root:
		#rootout = root(FPO, xn, args=(tof), method='broyden1', tol=1.49012e-08)
		#xpo = rootout.x
		#convergence = np.append(convergence, np.array([rootout.success], str))
		
		xpo = xpo.reshape(np.size(tof), 2)
		x2Dpo[indices[posn-1]] = xpo

	print "x2Dpo:"
	print x2Dpo
	
	convergence = convergence[1:]
	print "convergence:"
	print convergence
	
	xpo = np.array([ps2D2psxhat(x2D) for x2D in x2Dpo])
	print "xpo:"
	print xpo


	tofpo = np.array([timeofflight(x0) for x0 in xpo], float)
	T = np.zeros(np.max(group))
	
	for i in range(1,int(np.max(group))+1):
		#Get corresponding indices for the group i:
		indices = np.argwhere(group == i)
		#Reshape them in a 1D array:
		indices = indices.reshape(np.size(indices))
		T[i-1] = np.sum(tofpo[indices])
	
	print "Periods:"
	print T

	np.savetxt('data/group.dat', group)
	np.savetxt('data/itineraries.dat', scitineraries, fmt="%s")
	np.savetxt('data/position.dat', position)
	np.savetxt('data/tofpo.dat',tofpo)
	np.savetxt('data/xpo.dat', xpo)
	np.savetxt('data/periods.dat', T)

srange = np.linspace(np.min(sn), np.max(sn), 4000)
	
if plotpsandretmap:
	#Plotting:
	
	mpl.rcParams.update({'font.size': 22})

	
	lw = 3
	
	#Plot Poincare section:
	figure(1, figsize=(6, 6))
	
	xlabel('$e_1$', fontsize=36)
	ylabel('$e_2$', fontsize=36)
	
	plot(psectprojected[:,0],psectprojected[:,1], '.', ms=10)
	#plt.grid()
	plt.hold(True)
	
	plot(xintpsect, yintpsect, c='k', linewidth=lw)
	
	savefig('image/psectonslice.png', bbox_inches='tight', dpi=200)
	
	#Plot return maps:
	figure(2, figsize=(6, 6))
	
	xlabel('$s_n$', fontsize=36)
	ylabel('$s_{n+1}$', fontsize=36)
	
	plot(snsorted,snplus1sorted, '.', ms=10)
	#plt.grid()
	
	plt.hold(True)
	
	#sp1 = np.array([retmapm(1, sn) for sn in srange])
	
	plot(srange, srange, linewidth=lw)
	plot(snint, snplus1int, c='k', linewidth=lw)
	#plot(srange, sp1, c='k', linewidth=lw)
	
	savefig('image/retmaponslice.png', bbox_inches='tight', dpi=200)
	
	figure(3, figsize=(8,8))
	
	plt.subplot(231)	
	
	xlabel('$s_n$', fontsize=36)
	ylabel('$s_{n+2}$', fontsize=36)
	
	plt.hold(True)
	
	plot(srange, srange)
	
	sp2 = np.array([retmapm(2, sn) for sn in srange])
	
	plot(srange, sp2, c='k')
		
	plt.subplot(232)	
	
	xlabel('$s_n$', fontsize=36)
	ylabel('$s_{n+3}$', fontsize=36)
	
	plt.hold(True)
	
	plot(srange, srange)
	
	sp3 = np.array([retmapm(3, sn) for sn in srange])
	
	plot(srange, sp3, c='k')

	plt.subplot(233)
	
	xlabel('$s_n$', fontsize=36)
	ylabel('$s_{n+4}$', fontsize=36)
	
	plt.hold(True)
	
	plot(srange, srange)
	
	sp4 = np.array([retmapm(4, sn) for sn in srange])
	
	plot(srange, sp4, c='k')
	
	plt.subplot(234)
	
	xlabel('$s_n$', fontsize=36)
	ylabel('$s_{n+5}$', fontsize=36)
	
	plt.hold(True)
	
	plot(srange, srange)
	
	sp5 = np.array([retmapm(5, sn) for sn in srange])
	
	plot(srange, sp5, c='k')
	
	plt.subplot(235)
	
	xlabel('$s_n$', fontsize=36)
	ylabel('$s_{n+6}$', fontsize=36)
	
	plt.hold(True)
	
	plot(srange, srange)
	
	sp6 = np.array([retmapm(6, sn) for sn in srange])
	
	plot(srange, sp6, c='k')
	
	plt.subplot(236)
	
	xlabel('$s_n$', fontsize=36)
	ylabel('$s_{n+7}$', fontsize=36)
	
	plt.hold(True)
	
	plot(srange, srange)
	
	sp7 = np.array([retmapm(7, sn) for sn in srange])
	
	plot(srange, sp7, c='k')
	
	plt.tight_layout()
	plt.show()
