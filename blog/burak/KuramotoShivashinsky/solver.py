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

#Callable function version:
def integrate(x0, t, abserror=1.0e-8, relerror=1.0e-6):
	"""
	Takes the initial condition, parameters and the time interval
	returns the result as a series in time.
	ChaosBook Parameters [a=0.2, b=0.2, c=5.7] are assigned if parameters 
	are not specified otherwise
	"""
	xsol = odeint(KS.vfullssp, x0, t, atol = abserror, rtol = relerror)
	
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
	numsteps = 5000
	
	# Create the time samples for the output of the ODE solver:
	t = np.linspace(0, stoptime, numsteps+1)
	print "t"
	print t
	
	# Initial conditions
	
	x0 = [-0.0828769447817613, -0.08756936770734143, -0.5470748035590064, 
	0.1288931748248967, -0.07498656318075768, -0.04653195651665294, 0.129841924646638,
	0.3722174233509085, 0.01381517488153687, 0.04217832413769849, 0.05607306259312764,
	-0.0338005926473539, 0.009389167918579444, -0.001882241079403492, -0.007183362405766229,
	-0.009032458603270463, -0.0009662950988101977, -0.001712009441519379, 
	-0.001024726152019991, 0.001137241501704473, -0.0002493338657143959, 
	0.0001841065340444997, 0.0001605910827478338, 0.0001097200466209565, 
	3.412070953428844e-05, 3.161997837662267e-05, 1.042066241507e-05, 
	-1.993950459172108e-05, 4.555770323734273e-06, -5.736367915065645e-06]
	x0 = np.random.randn(30)
	
	# Call the ODE solver
	xsol = integrate(x0,t)

	# Write the solution into the file data/solution.dat
	# Use numpy.savetxt (practical, faster):
	tx = np.array([t], float).transpose() #Add time to the first collumn of tx
	tx = np.append(tx, xsol, axis=1) #Add solution to the following collumns of tx
	np.savetxt("data/sspsol.dat", tx) #Write tx in data/solution.dat

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
