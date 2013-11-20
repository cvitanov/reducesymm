#
# plotsliceandgorbit.py
#
"""
Plots the magic slice x'=[1,1,0,0] and the group orbit
"""

from __future__ import unicode_literals
from numpy import loadtxt
from pylab import figure, plot, xlabel, ylabel, grid, hold, legend, title, savefig
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np
import twomode

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)
        
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True

t, x1GS, y1GS, y2GS = loadtxt('data/gramschmidt.dat', unpack=True)

x1 = np.array([0.5,0.5,0.5,0.5], float)
x2 = np.array([0.1,0.1,1,1], float)
theta = np.arange(0, 2*np.pi, 0.01)
xg1	= np.zeros((np.size(theta,0),4))
xg2	= np.zeros((np.size(theta,0),4))

for i in range(0, np.size(theta,0), 1): 
	xg1[i,:] = np.dot(twomode.LieElement(theta[i]), x1)
	xg2[i,:] = np.dot(twomode.LieElement(theta[i]), x2)


fig = plt.figure()

ax = fig.gca(projection='3d')
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
ax.view_init(30,30)

x1 = [0,1,1,0]
x2 = [0,1,1,0]
y1 = [-2,-2,2,2]
verts = [zip(x1, x2, y1)]
poly = Poly3DCollection(verts, facecolor = '0.75', alpha=1)

ax.add_collection3d(poly)

ax.hold(True)

ax.plot(xg1[:,0], xg1[:,1], xg1[:,2], linewidth=1)
t1 = np.dot(twomode.generator(), xg1[0,:])/2
a1 = Arrow3D([xg1[0,0],xg1[0,0]+t1[0]],[xg1[0,1],xg1[0,1]+t1[1]],[xg1[0,2],xg1[0,2]+t1[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
ax.add_artist(a1)

ax.plot(xg2[:,0], xg2[:,1], xg2[:,2], linewidth=1, c='g')
#ax.plot(xg2[0,0], xg2[0,1], xg2[0,2], linewidth=1, c='g')
t2 = np.dot(twomode.generator(), xg2[0,:])/2
a2 = Arrow3D([xg2[0,0],xg2[0,0]+t2[0]],[xg2[0,1],xg2[0,1]+t2[1]],[xg2[0,2],xg2[0,2]+t2[2]], mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
ax.add_artist(a2)

ax.plot([xg1[0,0], xg2[0,0]], [xg1[0,1], xg2[0,1]], [xg1[0,2], xg2[0,2]], linestyle=' ', linewidth=1, marker='.', color='k')


ax.set_xlabel('$x_{1}$', fontsize=16)
ax.set_ylabel('$y_{1}$', fontsize=16)
ax.set_zlabel('$y_{1}$', fontsize=16)
savefig('image/sliceandgorbit.png', bbox_inches='tight', dpi=100)

plt.tight_layout()
plt.show()
