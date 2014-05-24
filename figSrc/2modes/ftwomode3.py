#
# invpolsolver.py --> ftwomode3.py
#
"""
2 mode figures in 4 different representations
"""

#from __future__ import unicode_literals

from scipy.integrate import odeint
import twomode
import numpy as np
import sspsolver3
import onslicesolver3
import invpolsolver3

import matplotlib as mpl
from pylab import figure, plot, xlabel, ylabel, grid, hold, legend, title, savefig
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from subprocess import call

#Load parameters:
p = np.loadtxt('data/parameters.dat')

x01 = [0.43997, 0, -0.38627, 0.07020] #reqv
x02 = [4.525719078826287434e-01, -2.791192295890519988e-18, 
5.092565177630036660e-02, 3.354280141917114627e-02] #rpo 1  
x03 = [0.0384074556708, 0.0, -1.90362452394, 0.0668631895808] #attractor

xphi01 = [0.43997, 0, -0.38627, 0.07020, 0] #reqv
xphi02 = [4.525719078826287434e-01, -2.791192295890519988e-18, 
5.092565177630036660e-02, 3.354280141917114627e-02, 0] #rpo 1 
xphi03 = [0.0384074556708, 0.0, -1.90362452394, 0.0668631895808, 0] #attractor

u01 = twomode.ssp2invpol(x01)
u02 = twomode.ssp2invpol(x02)
u03 = twomode.ssp2invpol(x03)

t1 = np.linspace(0,150,100000)
t2 = np.linspace(0, 10*3.641511999233241426e+00, 10000)
t3 = t1/2

xsol1 = sspsolver3.integrate(x01, p, t1)
xsol2 = sspsolver3.integrate(x02, p, t2)
xsol3 = sspsolver3.integrate(x03, p, t3)

xhatsol1 = onslicesolver3.integrate(xphi01, p, t1)
xhatsol2 = onslicesolver3.integrate(xphi02, p, t2)
xhatsol3 = onslicesolver3.integrate(xphi03, p, t3)

usol1 = invpolsolver3.integrate(u01, p, t1)
usol2 = invpolsolver3.integrate(u02, p, t2)
usol3 = invpolsolver3.integrate(u03, p, t3)

xtildesol1 = twomode.ssp2sspRed2(xsol1)
xtildesol2 = twomode.ssp2sspRed2(xsol2)
xtildesol3 = twomode.ssp2sspRed2(xsol3)

#Plot full state space:

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))

ax.plot(xsol3[:,0], xsol3[:,2], xsol3[:,3], linewidth=0.5, color='#3c5f96')
ax.hold(True)
ax.plot(xsol2[:,0], xsol2[:,2], xsol2[:,3], linewidth=0.8, color='#f7464a')
ax.plot(xsol1[:,0], xsol1[:,2], xsol1[:,3], linewidth=2.5, color='#33CC33')

ax.set_xlabel('\n $x_1$ \t', fontsize=32)
ax.set_ylabel('\n $x_2$ \t', fontsize=32)
ax.set_zlabel('$y_2$   ', fontsize=32)

Nticks = 3

xticks = np.linspace(-1.5, 1.5, Nticks)
ax.set_xticks(xticks)
ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=24); 

yticks = np.linspace(-1.5, 1.5, Nticks)
ax.set_yticks(yticks)
ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=24); 

zticks = np.linspace(-1.5, 1.5, Nticks)
ax.set_zticks(zticks)
ax.set_zticklabels(["$-1.5$ \t", " $0.0$", " $1.5$"], fontsize=24); 

ax.view_init(15,30)
savefig('2modes-ssp.pdf', bbox_inches='tight', dpi=100)
#savefig('twomode1.png', bbox_inches='tight', dpi=150)

#call(['convert', '-trim', "twomode1.png", "twomode1.png"])

#Plot reduced state space:
fig.clf()
ax = fig.gca(projection='3d')
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))

ax.plot(xhatsol3[:,0], xhatsol3[:,2], xhatsol3[:,3], linewidth=0.5, color='#3c5f96')
ax.hold(True)
ax.plot(xhatsol2[:,0], xhatsol2[:,2], xhatsol2[:,3], linewidth=3, color='#f7464a')
ax.plot(xhatsol1[:,0], xhatsol1[:,2], xhatsol1[:,3], linewidth=5, color='#33CC33')

ax.set_xlabel('\n $\hat{x}_1$ \t  ', fontsize=32)
ax.set_ylabel('\n $\hat{x}_2$ \t', fontsize=32)
ax.set_zlabel('$\hat{y}_2$   ', fontsize=32)

xticks = np.linspace(0, 2, Nticks)
ax.set_xticks(xticks) 
ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=24); 

yticks = np.linspace(-1.8, 0, Nticks)
ax.set_yticks(yticks)
ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=24); 

zticks = np.linspace(-0.3, 0.3, Nticks)
ax.set_zticks(zticks)
ax.set_zticklabels(["$%.1f$ " % ztik for ztik in zticks], fontsize=24); 

ax.view_init(15,30)
savefig('2modes-sspRed.pdf', bbox_inches='tight', dpi=100)
#savefig('twomode2.png', bbox_inches='tight', dpi=150)

#call(["pdfcrop", "twomode2.pdf", "twomode2.pdf"])
#call(['convert', '-trim', "twomode2.png", "twomode2.png"])

#Plot invariant polynomials:
fig.clf()
ax = fig.gca(projection='3d')
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))

ax.plot(usol3[:,0], usol3[:,1], usol3[:,2], linewidth=0.5, color='#3c5f96')
ax.hold(True)
ax.plot(usol2[:,0], usol2[:,1], usol2[:,2], linewidth=3, color='#f7464a')
ax.plot(usol1[:,0], usol1[:,1], usol1[:,2], linewidth=5, color='#33CC33')

ax.set_xlabel('\n $u$ \t  ', fontsize=32)
ax.set_ylabel('\n $v$ \t', fontsize=32)
ax.set_zlabel('$w$   ', fontsize=32)

xticks = np.linspace(0, 4, Nticks)
ax.set_xticks(xticks) 
ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=24); 

yticks = np.linspace(0, 4, Nticks)
ax.set_yticks(yticks)
ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=24); 

zticks = np.linspace(-10, 0, Nticks)
ax.set_zticks(zticks)
ax.set_zticklabels(["$%.1f$ " % ztik for ztik in zticks], fontsize=24); 

ax.view_init(15,30)
savefig('2modes-invpol.pdf', bbox_inches='tight', dpi=100)
#savefig('twomode2.png', bbox_inches='tight', dpi=150)

#call(["pdfcrop", "twomode2.pdf", "twomode2.pdf"])
#call(['convert', '-trim', "twomode2.png", "twomode2.png"])

#Plot 2nd mode slice:
fig.clf()
ax = fig.gca(projection='3d')
ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))

ax.plot(xtildesol3[:,0], xtildesol3[:,1], xtildesol3[:,2], linewidth=0.5, color='#3c5f96')
ax.hold(True)
ax.plot(xtildesol2[:,0], xtildesol2[:,1], xtildesol2[:,2], linewidth=3, color='#f7464a')
ax.plot(xtildesol1[:,0], xtildesol1[:,1], xtildesol1[:,2], linewidth=5, color='#33CC33')

ax.set_xlabel('\n \t   $\\tilde{x}_1$ ', fontsize=32)
ax.set_ylabel('\n $\\tilde{y}_1$   ', fontsize=32)
ax.set_zlabel(' $\\tilde{x}_2$\t', fontsize=32)

xticks = np.linspace(-1.5, 0.5, Nticks)
ax.set_xticks(xticks) 
ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=24); 

yticks = np.linspace(-1, -0.2, Nticks)
ax.set_yticks(yticks)
ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=24); 

zticks = np.linspace(0.2, 1.8, Nticks)
ax.set_zticks(zticks)
ax.set_zticklabels(["$%.1f$ " % ztik for ztik in zticks], fontsize=24); 

ax.view_init(30,-120)
savefig('2modes-sspRed2.pdf', bbox_inches='tight', dpi=100)
#savefig('twomode2.png', bbox_inches='tight', dpi=150)

#call(["pdfcrop", "twomode2.pdf", "twomode2.pdf"])
#call(['convert', '-trim', "twomode2.png", "twomode2.png"])

call(['bash', 'cropscript.sh'])

#plt.show()
