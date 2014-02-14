#
# diff.py
#
"""
Compute difference between reconstructed scaled time solution and the ETDRK4 solution 
"""

from __future__ import unicode_literals
import numpy as np

txscaledtime = np.loadtxt("data/solutionreconstructed.dat")

txETDRK4 = np.loadtxt("data/sspsolution.dat")

nofelements = 40000
tscaledn = txscaledtime[0:nofelements,0]

indices = np.zeros(nofelements, int) 

for i in range(nofelements):
    indices[i] = np.argmin(np.abs(txETDRK4[:,0] - tscaledn[i]))

diff = txscaledtime[0:nofelements,1:np.size(txscaledtime,1)] - txETDRK4[indices,1:np.size(txETDRK4,1)]
absdiff = np.linalg.norm(diff, axis=1)
normxST = np.linalg.norm(txscaledtime[range(nofelements), 1:np.size(txscaledtime,1)], axis=1)
reldiff = np.array([absdiff[i]/normxST[i] for i in range(nofelements)], float)

tETDRK4 = txETDRK4[indices,0]
tdiff = np.abs(tscaledn - tETDRK4)

#print diff
#print absdiff

from pylab import figure, grid, hold, legend, savefig, plot, xlabel, ylabel
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

plt.subplot(3, 1, 1)

plot(tscaledn, absdiff)
ylabel('$|x_{ETDRK4} - x_{ST}|$')

plt.subplot(3, 1, 2)

plot(tscaledn, reldiff)
ylabel('$|x_{ETDRK4} - x_{ST}|/|x_{ST}|$')

plt.subplot(3, 1, 3)

plot(tscaledn, normxST)
ylabel('$|t_{ETDRK4} - t_{ST}|$')

savefig("image/diff.png")

plt.show()
