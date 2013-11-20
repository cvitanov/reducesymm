#
# expandingdirections.py
#
"""
Computes the stability eigenvalues and eigenvectors for the origin and the
relative equilibrium and draws the directions towards which the flow will expand
in the symmetry reduced state space.
"""

from __future__ import unicode_literals
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

#Load the solution
t, x1hat, y1hat, x2hat, y2hat, phi = np.loadtxt('data/solutiononslice.dat', unpack=True)

#Load parameters:
p = np.loadtxt('data/parameters.dat')

x0 = [0,0,0,0]
A0 = twomode.StabilityMatrix(x0, p)
w0, v0 = np.linalg.eig(A0)

print 'Stability matrix at the origin'
print A0
print 'Stability eigenvalues at the origin:'
print w0
print 'Stability eigenvectors at the origin:'
print v0
vexp0 = w0[0]*v0[:,0]+w0[1]*v0[:,1]
vexp0 = vexp0/np.linalg.norm(vexp0)
print 'Expanding direction about the origin:'
print vexp0

#Stability matrix on slice is computed in /blog/burak/PorterKnobloch/Mathematica/PKrootsSet1.nb:

Ahatreqv = np.array([[-0.3871393533445726, -1.0881681068716151, -3.409733242959543,  0], 
				     [0, 0, 0, 0], 
				     [1.78403438176627, -1.910710833497562, 0.48510456780085, -1.17633621374323], 
				     [-0.16432158992199475, -9.632867689752878, 0.08816810687161492, -5.502034993628388]], float)

reqv = [0.4399655797367152, 0, -0.38626706847930564, 0.0702043939917171]

#Areqv = twomode.StabilityMatrix(reqv, p)
wreqv, vreqv = np.linalg.eig(Ahatreqv)

print 'Stability matrix at the relative equilibrium'
print Ahatreqv
print 'Stability eigenvalues at the relative equilibrium:'
print wreqv
print 'Stability eigenvectors at the relative equilibrium:'
print vreqv
treqv = np.dot(twomode.generator(), reqv)
treqvnormalized = treqv / np.linalg.norm(treqv)
print 'Normalized group tangent at the relative equilibrium:'
print treqvnormalized
vexpreqv = wreqv[0]*vreqv[:,0]+wreqv[1]*vreqv[:,1]
vexpreqv = vexpreqv/np.linalg.norm(vexpreqv)
print 'Expanding direction about the relative equilibrium:'
print vexpreqv

#Plot time evolution and unstable directions:
fig = plt.figure()

axfsize=16
lw = 0.6
ax = fig.gca(projection='3d')
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
ax.view_init(20,30)
ax.set_xlabel('$\hat{x}_1$', fontsize=axfsize)
ax.set_ylabel('$\hat{x}_2$', fontsize=axfsize)
ax.set_zlabel('$\hat{y}_2$', fontsize=axfsize)

ax.plot(x1hat, x2hat, y2hat, linewidth=0.5*lw)

vexp0 = np.real(vexp0)/1.5
ax.hold(True)
a1 = Arrow3D([0,0+vexp0[0]], [0,0+vexp0[2]], [0,0+vexp0[3]], mutation_scale=20, lw=2, arrowstyle="-|>", color="m")
ax.add_artist(a1)

vexpreqv = np.real(vexpreqv)/1.5
a2 = Arrow3D([reqv[0],reqv[0]+vexpreqv[0]], [reqv[2],reqv[2]+vexpreqv[2]], [reqv[3],reqv[3]+vexpreqv[3]], mutation_scale=20, lw=2, arrowstyle="-|>", color="r")
ax.add_artist(a2)

savefig('image/unstabledirections.png', bbox_inches='tight', dpi=100)

generatevtk = 0

if generatevtk:
	
	import subprocess
	subprocess.call("clear", shell=True)
	subprocess.call("python slicesol2vtk.py > vtk/solutiononslice.vtk")
	
	
	#process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	#output = process.communicate()[0]

else:
	plt.tight_layout()
	plt.show()
