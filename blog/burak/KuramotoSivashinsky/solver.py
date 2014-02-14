#
# invpolsolver.py
#
"""
Use odeint to solve differential equations defined by v in rossler.py
"""

from __future__ import unicode_literals
from scipy.integrate import odeint, ode
import KS
import numpy as np
import os
import oct2py


#Callable function version:
def integrate(x0, t, abserror=1.0e-8, relerror=1.0e-6):
	"""
	Takes the initial condition, parameters and the time interval
	returns the result as a series in time.
	ChaosBook Parameters [a=0.2, b=0.2, c=5.7] are assigned if parameters 
	are not specified otherwise
	"""
	
	mfilesdir = os.getcwd()
	oct2py.octave.addpath(mfilesdir)
	
	t, xsol = oct2py.octave.intwr5(t, x0)	
	
	#xsol = odeint(KS.vfullssp, x0, t, atol = abserror, rtol = relerror)
	
	#xsol = diffeq.pc4(KS.vfullssp, x0,t)
	
	#solver = odespy.ThetaRule(KS.vfullssp)
	#solver.set_initial_condition(x0)
	#xsol, t = solver.solve(t)
	
	#r = ode(KS.vfullssp).set_integrator('dopri5')#, method='bdf', with_jacobian=False)
	#r.set_initial_value(x0, t[0])
	#dt = t[1]-t[0]
	#xsol = np.zeros((np.size(t), np.size(x0)),float)
	#print "np.shape(xsol)"
	#print np.shape(xsol)
	#xsol[0,:] = x0 
	#i = 1
	#imax = np.size(t)
	
	#while r.successful() and i<imax:
		#r.integrate(r.t+dt)
		##print "i:"
		##print i
		##print "r.t"
		##print r.t
		#xsol[i,:] = r.y
		#i = i + 1
			
	return xsol

#If the solver is called as "main":
if __name__ == "__main__":
			
	stoptime = 50
	numsteps = 500
	
	# Create the time samples for the output of the ODE solver:
	t = np.linspace(0, stoptime, numsteps+1)
	#print "t"
	#print t
	
	# Initial conditions

	x0 = [1.081769255879405645e-01, 0.000000000000000000e+00, -1.130645021532460798e-01, 
	2.735234400271993951e-02, -2.300369007695817619e-02, 2.743365567776075153e-02,
	4.242109469805525057e-01, -3.221375201839944691e-02, 3.517350195620121411e-01,
	4.196323928512094570e-01, 7.405822221667555938e-02, -4.911698645198029345e-01,
	-2.619084037151255262e-01, 8.869647954573157966e-03, 2.667068425090810685e-02,
	-1.245857190912121742e-01, 1.848625450932936676e-01, -1.886910780372257068e-01,
	-4.364329592632312099e-02, -8.603322827952401136e-03, -4.893648586116418342e-02,
	-4.227361593906937137e-02, -5.743046880645331920e-02, 6.141174799787345318e-02, 
	3.556320072237982056e-03, -2.647610106987533310e-02, -3.295731184608969265e-03, 
	-1.760410218329051119e-02, -1.449156681256755577e-02, 1.551927466950007474e-02]
	
	#x0 = np.random.randn(30)
	
	# Call the ODE solver
	#xsol = integrate(x0,t)
	xsol = integrate(x0,[0, stoptime])

	# Write the solution into the file data/solution.dat
	# Use numpy.savetxt (practical, faster):
	#tx = np.array([t], float).transpose() #Add time to the first collumn of tx
	#tx = np.append(tx, xsol, axis=1) #Add solution to the following collumns of tx
	#np.savetxt("data/sspsol.dat", tx) #Write tx in data/solution.dat

	# Plot the solution:
	#Import plotting modules:
	from pylab import figure, grid, hold, legend, savefig
	import matplotlib as mpl
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt 

	mpl.rcParams['text.usetex']=True
	mpl.rcParams['text.latex.unicode']=True

	fig = plt.figure() #Create a figure instance
	ax = fig.gca(projection = "3d") #Create an axis instance
	
	ax.plot(xsol[:,0], xsol[:,1], xsol[:,2], lw=0.5)
	ax.set_xlabel('$x$', fontsize=18)
	ax.set_ylabel('$y$', fontsize=18)
	ax.set_zlabel('$z$', fontsize=18)
	savefig("image/sspsolodeint.png", bbox_inches='tight', dpi=100)

	plt.tight_layout()
	plt.show()
