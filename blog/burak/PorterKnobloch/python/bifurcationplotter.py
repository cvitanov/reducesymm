#
# psectgsplotter.py
#
"""
Plot the Poincar\'e section
"""
from __future__ import unicode_literals
from numpy import loadtxt
from pylab import figure, plot, xlabel, ylabel, xlim, grid, hold, legend, title, savefig
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

x = loadtxt('data/bifurcation.dat', unpack=True)
x=x.transpose()

figure(1, figsize=(5, 5))

plt.subplot(2,1,1)

plt.plot(x[:,0], x[:,1], '.', ms=0.1)
plt.ylabel('$\hat{x}_1,GS$')
xlim([np.min(x[:,0])-np.min(x[:,0])/100, np.max(x[:,0])+np.max(x[:,0])/100])
plt.grid()

plt.subplot(2,1,2)

plt.plot(x[:,0], x[:,3], '.', ms=0.1)
xlim([np.min(x[:,0])-np.min(x[:,0])/100, np.max(x[:,0])+np.max(x[:,0])/100])
plt.grid()

plt.ylabel('$\hat{y}_2,GS$')
plt.xlabel('$a_1$')

savefig('image/bifurcation.png', bbox_inches='tight', dpi=200)
plt.tight_layout()
plt.show()
