import numpy as np
import matplotlib as mpl
from pylab import plot, xlabel, ylabel, show, savefig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from subprocess import call

EscRateZeta = np.loadtxt('data/escratedynzeta.dat')
EscRateSpecDet = np.loadtxt('data/escratespecdet.dat')

EscRateZeta = -EscRateZeta
EscRateSpecDet = -EscRateSpecDet

gammaInfZeta = EscRateZeta[-1].copy()
gammaInfSpecDet = EscRateSpecDet[-1].copy()

ConvergenceZeta = np.log(np.abs(EscRateZeta[0:-1] - gammaInfZeta))
ConvergenceSpecDet = np.log(np.abs(EscRateSpecDet[0:-1] - gammaInfSpecDet))

fig=plt.figure(1, figsize=(8,6))
nrange = range(1,10)
plot(nrange, ConvergenceZeta, lw=2)
plt.hold('True')
plot(nrange, ConvergenceSpecDet, 'r', lw=2)
xlabel('$N$', fontsize=24)
ylabel('$\ln|\gamma_n - \gamma_{\infty}|$', fontsize=24)

Nticks = 9
ax=fig.gca()
ax.set_xticks(nrange)
ax.set_xticklabels(["$%1i$" % xtik for xtik in nrange], fontsize=16);
yticks = np.linspace(0, -16, Nticks)
ax.set_yticks(yticks)
ax.set_yticklabels(["$%2i$" % ytik for ytik in yticks], fontsize=16);

savefig('GammaConv.pdf', bbox_inches='tight', dpi=100)
call(["pdfcrop", "GammaConv.pdf", "GammaConv.pdf"], shell=True)

plt.show()
