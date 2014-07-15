#
# invpolsolver.py
#
"""
2 mode figures in 4 different representations
"""
from scipy.integrate import odeint
import twomode
import numpy as np
import sspsolver
import onslicesolver
import invpolsolver

import matplotlib as mpl
from pylab import figure, plot, xlabel, ylabel, grid, hold, legend, title, savefig, imshow
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from subprocess import call

#Load parameters:
p = np.loadtxt('data/parameters.dat')

x01 = [0.43997, 0, -0.38627, 0.07020] #reqv
x02 = [4.525719078826287434e-01, -2.791192295890519988e-18, 
5.092565177630036660e-02, 3.354280141917114627e-02] #rpo1  
x02 = [0.43998243, 0.0, -1.68577368, 0.06639063]# rpo01

x03 = [0.0384074556708, 0.0, -1.90362452394, 0.0668631895808] #attractor

xphi01 = [0.43997, 0, -0.38627, 0.07020, 0] #reqv
xphi02 = [4.525719078826287434e-01, -2.791192295890519988e-18, 
5.092565177630036660e-02, 3.354280141917114627e-02, 0] #rpo 1 
xphi02 = [0.43998243, 0.0, -1.68577368, 0.06639063, 0] #rpo 01 
xphi03 = [0.0384074556708, 0.0, -1.90362452394, 0.0668631895808, 0] #attractor

u01 = twomode.ssp2invpol(x01)
u02 = twomode.ssp2invpol(x02)
u03 = twomode.ssp2invpol(x03)

tfreqv = 50
tfrpo = 10*3.641511999233241426e+00
tfrpo = 2* 7.34594127
tfergo = tfreqv
treqv = np.linspace(0,tfreqv,5000)
trpo = np.linspace(0, tfrpo, 10000)
tergo = np.linspace(0, tfergo, 5000)

xsolreqv = sspsolver.integrate(x01, p, treqv)
xsolrpo = sspsolver.integrate(x02, p, trpo)
xsolergo = sspsolver.integrate(x03, p, tergo)

xhatsolreqv = onslicesolver.integrate(xphi01, p, treqv)
xhatsolrpo = onslicesolver.integrate(xphi02, p, trpo)
xhatsolergo = onslicesolver.integrate(xphi03, p, tergo)


def inverseFourier(xsol):
    """
    Generates configuration space solution u(x) = F^-1{z[k]} from the Fourier 
    modes
    """
    zsol = np.array([[xsol[i,0] + 1j*xsol[i,1], xsol[i,2] + 1j*xsol[i,3]] for i 
                    in range(np.size(xsol, 0))],complex)
    
    x = np.arange(-np.pi, np.pi, 0.05)
    exp1jx = np.exp(1j*x) 
    exp2jx = np.exp(2j*x) 
    usol = np.array([[2*np.real(zsol[i,0]*exp1jx) + 
                      2*np.real(zsol[i,1]*exp2jx)]
                      for i in range(np.size(zsol,0))], float)
    usol = usol.reshape((np.size(usol,0), np.size(usol, 2)))
    return usol

#Compute inverse Fourier transforms of the equivariant solutions
usolreqv = inverseFourier(xsolreqv)
usolrpo = inverseFourier(xsolrpo)
usolergo = inverseFourier(xsolergo)
#Compute inverse Fourier transforms of the invariant solutions
uhatsolreqv = inverseFourier(xhatsolreqv[:,0:4])
uhatsolrpo = inverseFourier(xhatsolrpo[:,0:4])
uhatsolergo = inverseFourier(xhatsolergo[:,0:4])

fig = plt.figure(figsize=(3, 6)) #Create a figure instance

x = np.arange(-np.size(usolreqv, 1)/2, np.size(usolreqv, 1)/2)
y = treqv[range(0,np.size(treqv),10)]

im = plt.pcolormesh(x, y, usolreqv[range(0,np.size(treqv),10),:], shading='gouraud')

plt.axis([x.min(), x.max(), y.min(), y.max()])

plt.xticks([x[0],x[np.floor(np.size(x)/2)],x[-1]], ('$-\pi / 2$', '$0$', '$\pi / 2$'), fontsize=32)
plt.yticks([0,tfreqv/2.0,tfreqv], 
('$'+str(0)+'$', '$'+str(tfreqv/2)[0:2]+'$', '$'+str(tfreqv)[0:2]+'$'), fontsize=32)
plt.xlabel('$x$', fontsize=40)
plt.ylabel('$t$', fontsize=40)

savefig('2modes-conf-reqv.png', bbox_inches='tight', dpi=100)
call(['convert', '-trim', '2modes-conf-reqv.png', '2modes-conf-reqv.png'])

fig.clf()

x = np.arange(-np.size(uhatsolreqv, 1)/2, np.size(uhatsolreqv, 1)/2)
y = treqv[range(0,np.size(treqv),10)]

im = plt.pcolormesh(x, y, uhatsolreqv[range(0,np.size(treqv),10),:], shading='gouraud')

plt.axis([x.min(), x.max(), y.min(), y.max()])

plt.xticks([x[0],x[np.floor(np.size(x)/2)],x[-1]], ('$-\pi / 2$', '$0$', '$\pi / 2$'), fontsize=32)
plt.yticks([0,tfreqv/2.0,tfreqv], 
('$'+str(0)+'$', '$'+str(tfreqv/2)[0:2]+'$', '$'+str(tfreqv)[0:2]+'$'), fontsize=32)
plt.xlabel('$x$', fontsize=40)
plt.ylabel('$t$', fontsize=40)

savefig('2modes-confred-reqv.png', bbox_inches='tight', dpi=100)
call(['convert', '-trim', '2modes-confred-reqv.png', '2modes-confred-reqv.png'])

fig.clf()

x = np.arange(-np.size(usolrpo, 1)/2, np.size(usolrpo, 1)/2)
y = trpo[range(0,np.size(trpo),10)]

im = plt.pcolormesh(x, y, usolrpo[range(0,np.size(trpo),10),:], shading='gouraud')

plt.axis([x.min(), x.max(), y.min(), y.max()])

plt.xticks([x[0],x[np.floor(np.size(x)/2)],x[-1]], ('$-\pi / 2$', '$0$', '$\pi / 2$'), fontsize=32)
plt.yticks([0,tfrpo/2,tfrpo], 
('$'+str(0)+'$', '$'+str(tfrpo/2.0)[0:3]+'$', '$'+str(tfrpo)[0:4]+'$'), fontsize=32)
plt.xlabel('$x$', fontsize=40)
plt.ylabel('$t$', fontsize=40)

savefig('2modes-conf-rpo.png', bbox_inches='tight', dpi=100)
call(['convert', '-trim', '2modes-conf-rpo.png', '2modes-conf-rpo.png'])

fig.clf()

x = np.arange(-np.size(uhatsolrpo, 1)/2, np.size(uhatsolrpo, 1)/2)
y = trpo[range(0,np.size(trpo),10)]

im = plt.pcolormesh(x, y, uhatsolrpo[range(0,np.size(trpo),10),:], shading='gouraud')

plt.axis([x.min(), x.max(), y.min(), y.max()])

plt.xticks([x[0],x[np.floor(np.size(x)/2)],x[-1]], ('$-\pi / 2$', '$0$', '$\pi / 2$'), fontsize=32)
plt.yticks([0,tfrpo/2,tfrpo], 
('$'+str(0)+'$', '$'+str(tfrpo/2.0)[0:3]+'$', '$'+str(tfrpo)[0:4]+'$'), fontsize=32)
plt.xlabel('$x$', fontsize=40)
plt.ylabel('$t$', fontsize=40)

savefig('2modes-confred-rpo.png', bbox_inches='tight', dpi=100)
call(['convert', '-trim', '2modes-confred-rpo.png', '2modes-confred-rpo.png'])

fig.clf()

x = np.arange(-np.size(usolergo, 1)/2, np.size(usolergo, 1)/2)
y = tergo[range(0,np.size(tergo),10)]

im = plt.pcolormesh(x, y, usolergo[range(0,np.size(tergo),10),:], shading='gouraud')

plt.axis([x.min(), x.max(), y.min(), y.max()])

plt.xticks([x[0],x[np.floor(np.size(x)/2)],x[-1]], ('$-\pi / 2$', '$0$', '$\pi / 2$'), fontsize=32)
plt.yticks([0,tfergo/2.0,tfergo], 
('$'+str(0)+'$', '$'+str(tfergo/2.0)[0:2]+'$', '$'+str(tfergo)+'$'), fontsize=32)
plt.xlabel('$x$', fontsize=40)
plt.ylabel('$t$', fontsize=40)

savefig('2modes-conf-ergodic.png', bbox_inches='tight', dpi=100)
call(['convert', '-trim', '2modes-conf-ergodic.png', '2modes-conf-ergodic.png'])

fig.clf()

x = np.arange(-np.size(uhatsolergo, 1)/2, np.size(uhatsolergo, 1)/2)
y = tergo[range(0,np.size(tergo),10)]

im = plt.pcolormesh(x, y, uhatsolergo[range(0,np.size(tergo),10),:], shading='gouraud')

plt.axis([x.min(), x.max(), y.min(), y.max()])

plt.xticks([x[0],x[np.floor(np.size(x)/2)],x[-1]], ('$-\pi / 2$', '$0$', '$\pi / 2$'), fontsize=32)
plt.yticks([0,tfergo/2.0,tfergo], 
('$'+str(0)+'$', '$'+str(tfergo/2.0)[0:2]+'$', '$'+str(tfergo)+'$'), fontsize=32)
plt.xlabel('$x$', fontsize=40)
plt.ylabel('$t$', fontsize=40)

savefig('2modes-confred-ergodic.png', bbox_inches='tight', dpi=100)
call(['convert', '-trim', '2modes-confred-ergodic.png', '2modes-confred-ergodic.png'])
