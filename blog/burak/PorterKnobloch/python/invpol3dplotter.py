#
# plotter3d.py
#
"""
Plot the solution generated by invpolsolver.py
"""

from numpy import loadtxt
from pylab import figure, plot, xlabel, ylabel, grid, hold, legend, title, savefig
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

t, u, v, w, q = loadtxt('data/invpolsolution.dat', unpack=True)

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(u, v, w, linewidth=0.3)
#ax.grid(False)
ax.set_xlabel('u')
ax.set_ylabel('v')
ax.set_zlabel('w')
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
ax.view_init(15,30)
savefig('image/uvw.png', bbox_inches='tight', dpi=100)

reqv = 0

fig2 = plt.figure()
plt.plot(t[95000:100000],q[95000:100000])
plt.xlabel('Time')
plt.ylabel('u')
savefig('image/qt.png', bbox_inches='tight', dpi=100)

if reqv:
	
	ax.hold(True)
	ax.plot(u[0:10000], v[0:10000], w[0:10000], linewidth=0.3, c='r')
	savefig('image/uvwreqv.png', bbox_inches='tight', dpi=150)

plt.show()
