#
# periodicorbits.py
#
import numpy as np

import onslicesolver
pars = np.loadtxt('data/parameters.dat')

xpo = np.loadtxt('data/xrpo.dat')
itineraries = np.loadtxt('data/itineraries.dat', dtype="str")
periods = np.loadtxt('data/periods.dat')
position = np.loadtxt('data/position.dat')
group = np.loadtxt('data/group.dat')
tofpo = np.loadtxt('data/tofrpo.dat')
print "group:"
print group
#Plotting modules:
from pylab import figure, grid, hold, legend, savefig
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

#fig = plt.figure() #Create a figure instance
#ax = fig.gca(projection = "3d") #Create an axis instance
#hold(False)

for i in range(1,int(np.max(group))+1):
	
	fig = plt.figure() #Create a figure instance
	ax = fig.gca(projection = "3d") #Create an axis instance
		
	#indices corresponding to the members of ith periodic orbit:
	gindices = np.argwhere(group == i)
	gindices = gindices.reshape(np.size(gindices))
	print "gindices:"
	print gindices
	#Position of first element of the cycle:
	p1 = np.argwhere(position[gindices] == 1)
	p1 = p1.reshape(np.size(p1))
	print "p1:"
	print p1
	#Position of the first element:
	i1 = gindices[p1]
	i1 = i1.reshape(np.size(i1))
	print "i1:"
	print i1
	x0 = xpo[i1,:]
	x0 = x0.reshape(np.size(x0))
	print "x0:"
	print x0
	itinerary = itineraries[i1]
	print itinerary[0]
	
	#Period of the ith cycle:
	T = periods[i-1]

	stoptime = T
	numpoints = int(stoptime / 0.01)
	
	# Create the time samples for the output of the ODE solver:
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	
	# Call the ODE solver
	xhatphi0 = np.append(x0, np.array([0], float))
	xhatsol = onslicesolver.integrate(xhatphi0, pars, t, abserror=1.0e-14, relerror=1.0e-12)
	
	ax.plot(xhatsol[:,0], xhatsol[:,2], xhatsol[:,3], lw=0.5)
	ax.set_xlabel('$x$', fontsize=18)
	ax.set_ylabel('$y$', fontsize=18)
	ax.set_zlabel('$z$', fontsize=18)
	#ax.view_init(30,-110)
	ax.grid(b='off')
	
	
	fname = 'image/' + itinerary[0] + '.png'
	savefig(fname, bbox_inches='tight', dpi=100)	

	plt.show()	
	raw_input("Press Enter to continue...")
