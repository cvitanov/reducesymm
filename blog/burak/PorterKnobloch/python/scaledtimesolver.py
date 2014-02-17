#
# invpolsolver.py
#
"""
Use odeint to solve differential equations defined by v in rossler.py
"""

from __future__ import unicode_literals
from scipy.integrate import odeint
import twomode
import numpy as np

#Callable function version:
def integrate(x0, t, p, abserror=1.0e-8, relerror=1.0e-6):
	"""
	Takes the initial condition, parameters and the time interval
	returns the result as a series in time.
	ChaosBook Parameters [a=0.2, b=0.2, c=5.7] are assigned if parameters 
	are not specified otherwise
	"""
	xsol = odeint(twomode.vscaledtime, x0, t, args=(p,), atol = abserror, rtol = relerror)
	
	return xsol

#If the solver is called as "main":
if __name__ == "__main__":
	
	p = np.loadtxt('data/parameters.dat')
			
	stoptime = 100000
	#stoptime =  35.06091404
	numpoints = 10000000
	#numpoints = int(stoptime / 0.01)
	
	# Create the time samples for the output of the ODE solver:
	#tau = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	tau = np.linspace(0,stoptime, numpoints)
	
	# Pack up the initial conditions (picked on the strange attractor):
	x0 = [0.43997, 0, -0.38627, 0.07020]
	#x0 =  [10.1066516,   -5.83507802,   0.54549834]
	
	# Call the ODE solver
	xsol = integrate(x0,tau,p)
	#print xsol
	
	# Write the solution into the file data/solution.dat
	# Use numpy.savetxt (practical, faster):
	taux = np.array([tau], float).transpose() #Add time to the first collumn of tx
	taux = np.append(taux, xsol, axis=1) #Add solution to the following collumns of tx
	#np.savetxt("data/solutionscaledtime.dat", taux) #Write tx in data/solution.dat
	
	#Use file module (more general, enables custom formatting):
	#f = open("data/solution.dat", "w")
	
	#for i in range(len(t)):
		
		#f.write("%5.7f " % float(t[i])) # Write time
		#for j in range(np.size(xsol, 1)):
			#f.write("\t %5.7f" % float(xsol[i][j])) #Write coordinates delimited by \t characters
		#f.write("\n") #New line

	#f.close()

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
	
	ax.plot(xsol[:,0], xsol[:,2], xsol[:,3], lw=0.3)
	ax.set_xlabel('$\hat{x}_1$', fontsize=18)
	ax.set_ylabel('$\hat{x}_2$', fontsize=18)
	ax.set_zlabel('$\hat{y}_2$', fontsize=18)
	ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
	ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
	ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
	ax.view_init(25,35)
	savefig("image/scaledtimesol.eps", bbox_inches='tight', dpi=100)

	plt.tight_layout()
	plt.show()
