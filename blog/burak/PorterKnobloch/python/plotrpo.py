#
# periodicorbits.py
#
import numpy as np

import onslicesolver
import sspsolver

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
from pylab import figure, grid, hold, legend, savefig, imshow, colorbar
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
	print "Period ="
	print T

	stoptime = T
	numpoints = int(stoptime / 0.01)
	
	# Create the time samples for the output of the ODE solver:
	t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	
	# Call the ODE solver
	xhatphi0 = np.append(x0, np.array([0], float))
	xhatsol = onslicesolver.integrate(xhatphi0, pars, t, abserror=1.0e-12, relerror=1.0e-10)
	
	ax.plot(xhatsol[:,0], xhatsol[:,2], xhatsol[:,3], lw=0.5)
	ax.set_xlabel('$\hat{x}_1$', fontsize=18)
	ax.set_ylabel('$\hat{x}_2$', fontsize=18)
	ax.set_zlabel('$\hat{y}_2$', fontsize=18)
	ax.view_init(25,35)
	ax.grid(b='off')
	
	
	fname = 'image/' + itinerary[0] + 'slice.png'
	savefig(fname, bbox_inches='tight', dpi=100)	
	
	plt.close()
	
	fig = plt.figure() #Create a figure instance
	ax = fig.gca(projection = "3d") #Create an axis instance
	
	stoptime = 2*T
	numpoints = int(stoptime / 0.01)
	
	# Create the time samples for the output of the ODE solver:
	#t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	
	t = np.linspace(0,stoptime,1001)
	#print "np.size(t,0)"
	#print np.size(t,0)
	xsol = sspsolver.integrate(x0, pars, t, abserror=1.0e-8, relerror=1.0e-6)
	
	ax.plot(xsol[:,0], xsol[:,2], xsol[:,3], lw=0.5)
	ax.set_xlabel('$x_1$', fontsize=18)
	ax.set_ylabel('$y_1$', fontsize=18)
	ax.set_zlabel('$x_2$', fontsize=18)
	ax.view_init(25,35)
	ax.grid(b='off')
	
	fname = 'image/' + itinerary[0] + 'ssp.png'
	savefig(fname, bbox_inches='tight', dpi=100)
	
	plt.close()
	
	#Complex Fourier coefficients:
	a1 = np.array([xsol[i,0] + xsol[i,1]*1j for i in range(np.size(xsol,0))], complex)
	am1 = np.conjugate(a1)
	a2 = np.array([xsol[i,2] + xsol[i,3]*1j for i in range(np.size(xsol,0))], complex)
	am2 = np.conjugate(a2)
	
	#From Fourier coefficients to the real space definition
	
	x = np.linspace(-np.pi,np.pi,501)#np.size(a1,0))
	u = np.zeros((np.size(xsol,0), np.size(x,0)), complex)
	
	for k in range(np.size(u, 0)):
		for l in range(np.size(u, 1)):
			u[k,l] = a1[k]*np.exp(1j*x[l]) + am1[k]*np.exp(-1j*x[l]) + a2[k]*np.exp(2j*x[l]) + am2[k]*np.exp(-2j*x[l])
	
	u=np.real(u)
	#Plot with imshow:
	#http://stackoverflow.com/questions/11775354/how-can-i-display-a-np-array-with-pylab-imshow
	
	fig = plt.figure() #Create a figure instance
	#ax = fig.gca() #Create an axis instance
	
	#im = ax.imshow(u, origin='lower')#, cmap='hot')
	
	im = imshow(u, origin='lower')
	colorbar(im, orientation='vertical')
	
	plt.xticks([0,250,501], ('$-\pi$', '$0$', '$-\pi$'))
	plt.xlabel('$space$')
	plt.yticks([0,500,1001], ('$0$', '$T$', '$2T$'))
	plt.ylabel('$time$')
	
	#xtickx = np.linspace(np.min(x), np.max(x), 9)
	#ttickt = np.linspace(np.min(t), np.max(t), 9)
	#print "xtickx"
	#print xtickx
	#ax.set_xticklabels(["%.2f" % xtick for xtick in xtickx])
	#ax.set_yticklabels(["%.2f" % ytick for ytick in ttickt])
	
	fname = 'image/' + itinerary[0] + 'conf.png'
	savefig(fname, bbox_inches='tight', dpi=200)

	#plt.show()	
	#raw_input("Press Enter to continue...")
	
	plt.close()
	

	#plt.show()	
	#raw_input("Press Enter to continue...")
