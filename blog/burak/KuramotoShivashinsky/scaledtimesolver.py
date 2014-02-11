#
# invpolsolver.py
#
"""
Use odeint to solve differential equations defined by v in rossler.py
"""

from __future__ import unicode_literals
from scipy.integrate import odeint
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
	xsol = odeint(KS.vscaledtime, x0, t, atol = abserror, rtol = relerror)
	
	return xsol

#If the solver is called as "main":
if __name__ == "__main__":
			
	stoptime = 100
	numpoints = 1000
	
	# Create the time samples for the output of the ODE solver:
	#tau = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	tau = np.linspace(0,stoptime, numpoints)
	
	# Initial conditions:
	#x0 = [1.205694079649668737e-01, 0.000000000000000000e+00, 1.587973000268175505e-01,
	#5.391546241197969769e-01, -2.713877446516336014e-02, -8.397436844334366102e-02,
	#-1.699516623357091305e-01, -3.556981418439826848e-01, 4.197035769714546616e-02,
	#1.443465132401413990e-02, 4.255859995423756648e-02, 4.975373334555246996e-02,
	#-8.839003381658995864e-03, -3.683927917766105335e-03, -8.982723411097105759e-03, 
	#-7.245459543196561616e-03, 1.965326008508331938e-03, 4.686391060090486227e-05, 
	#1.372928645112300641e-03, 6.770885136101844996e-04, -3.063006329784660471e-04, 
	#4.735519725387917476e-05, -1.874944200907897548e-04, -5.170906036126518974e-05,
	#4.292240638122765982e-05, -1.793635643340408189e-05, 2.239418180365538495e-05,
	#2.162098489187738860e-06, -6.331242075178413828e-06, 3.684607697585050442e-06]
	
	x0 = np.random.randn(62)
	x0[0] = np.abs(x0[0])
	x0[1] = 0
	
	# Call the ODE solver
	xsol = integrate(x0,tau)
	#print xsol
	
	# Write the solution into the file data/solution.dat
	# Use numpy.savetxt (practical, faster):
	taux = np.array([tau], float).transpose() #Add time to the first collumn of tx
	taux = np.append(taux, xsol, axis=1) #Add solution to the following collumns of tx
	np.savetxt("data/solutionscaledtime.dat", taux) #Write tx in data/solution.dat
	
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
	
	ax.plot(xsol[:,0], xsol[:,2], xsol[:,3], lw=0.5)
	ax.set_xlabel('$\hat{x}_1$', fontsize=18)
	ax.set_ylabel('$\hat{x}_2$', fontsize=18)
	ax.set_zlabel('$\hat{y}_2$', fontsize=18)
	ax.view_init(20,120)
	savefig("image/scaledtimesol.png", bbox_inches='tight', dpi=100)

	plt.tight_layout()
	plt.show()
