#
# psectgsplotter.py
#
"""
Plot the Poincar\'e section
"""
from __future__ import unicode_literals
from numpy import loadtxt
from pylab import figure, plot, xlabel, ylabel, grid, hold, legend, title, savefig
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

x = loadtxt('data/psectgs.dat', unpack=True)
x=x.transpose()

figure(1, figsize=(6, 4.5))
xlabel('$\hat{x}_1,GS$')
ylabel('$\hat{y}_2,GS$')

plot(x[:,1], x[:,3], '.', ms=2)
plt.grid()

savefig('image/psectgs.png', bbox_inches='tight', dpi=150)

plt.tight_layout()
plt.show()
