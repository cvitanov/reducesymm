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
			
	stoptime = 5000
	numpoints = 50000
	
	# Create the time samples for the output of the ODE solver:
	#tau = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]
	tau = np.linspace(0,stoptime, numpoints)
	
	# Initial conditions:
	x0 = [1.081769255879405645e-01, 0.000000000000000000e+00, -1.130645021532460798e-01, 
	2.735234400271993951e-02, -2.300369007695817619e-02, 2.743365567776075153e-02,
	4.242109469805525057e-01, -3.221375201839944691e-02, 3.517350195620121411e-01,
	4.196323928512094570e-01, 7.405822221667555938e-02, -4.911698645198029345e-01,
	-2.619084037151255262e-01, 8.869647954573157966e-03, 2.667068425090810685e-02,
	-1.245857190912121742e-01, 1.848625450932936676e-01, -1.886910780372257068e-01,
	-4.364329592632312099e-02, -8.603322827952401136e-03, -4.893648586116418342e-02,
	-4.227361593906937137e-02, -5.743046880645331920e-02, 6.141174799787345318e-02, 
	3.556320072237982056e-03, -2.647610106987533310e-02, -3.295731184608969265e-03, 
	-1.760410218329051119e-02, -1.449156681256755577e-02, 1.551927466950007474e-02, 0, 0]
	
	#x0 = np.random.randn(30)
	#x0[0] = np.abs(x0[0])
	#x0[1] = 0
	
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
