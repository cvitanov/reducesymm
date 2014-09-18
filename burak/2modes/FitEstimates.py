import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl
from pylab import plot, xlabel, ylabel, show, savefig
import matplotlib.pyplot as plt

EscapeRate = np.loadtxt('data/EscapeRate.dat')
AveragePeriod = np.loadtxt('data/AveragePeriod.dat')
AveragePhaseSpeed = np.loadtxt('data/AveragePhaseSpeed.dat')
DiffusionCoefficient = np.loadtxt('data/DiffusionCoefficient.dat')
LyapunovExponent = np.loadtxt('data/AverageLyapunovExponent.dat')

Nexpansion = 10

def Exponential(x, a,b,c):
    return a*np.exp(-b*x) + c

xdataeven = np.array([2], int)
i = 4
while i <= Nexpansion:
    xdataeven = np.append(xdataeven, i) 
    i += 2
    
xdataodd = np.array([1], int)
i = 3
while i <= Nexpansion:
    xdataodd = np.append(xdataodd, i) 
    i += 2

f = open("tex/DynamicalAverages.tex", "w")

f.write("\\begin{table}\n")
f.write("\t\\begin{tabular}{c|c|c|c|c|c}\n")

f.write("\t $N$ & $\gamma$ & $\langle T \\rangle$ & $\lambda$ & $\langle \
\dot{\phi} \\rangle$ & $D$ \\\\ \n")
f.write("\t\\hline\n")

for i in range(len(EscapeRate)):
    f.write("\t%s & " % str(i+1))
    f.write("%5.9f & " % float(EscapeRate[i]))
    f.write("%5.7f & " % float(AveragePeriod[i]))
    f.write("%5.8f & " % float(LyapunovExponent[i]))
    f.write("%5.7f & " % float(AveragePhaseSpeed[i]))
    f.write("%5.6f \\\\ \n " % float(DiffusionCoefficient[i]))
    
f.write("\t\\end{tabular}\n")
f.write("\t\\caption{Cyle expansion estimates of the escape rate $\gamma$, \
average cycle period $\langle T \\rangle$, Lyapunov exponent $\lambda$, \
average phase velocity $\langle \dot{\phi} \\rangle$ and the diffusion \
coefficient $D$ with respect to the expansion order $N$ .}\n")
f.write("\t\\label{t-DynamicalAverages}\n")
f.write("\\end{table}")

f.close()

#poptPeriodeven, pcovPeriododd = curve_fit(Exponential, xdataeven, AveragePeriod[xdataeven-1])
#poptPeriododd, pcovPeriododd = curve_fit(Exponential, xdataodd, AveragePeriod[xdataodd-1])
#aPeriododd, bPeriododd, cPeriododd = poptPeriododd
#aPeriodeven, bPeriodeven, cPeriodeven = poptPeriodeven

#plt.figure(1)
#plot(range(1,Nexpansion+1), AveragePeriod)
#plt.hold('True')
#xRange = np.arange(1,Nexpansion+0.01, 0.01)
#plot(xRange, [Exponential(xRange[i], aPeriododd, bPeriododd, cPeriododd) \
#              for i in range(len(xRange))])
#plot(xRange, [Exponential(xRange[i], aPeriodeven, bPeriodeven, cPeriodeven) \
#              for i in range(len(xRange))])

#poptPhase, pcovPhase = curve_fit(Exponential, xdata, AveragePhase[discard:len(AveragePhase)])

#aPhase, bPhase, cPhase = poptPhase
#aPeriod, bPeriod, cPeriod = poptPeriod
#plt.figure(1)
#xinterp = np.arange(0.0,float(len(AveragePhase)-discard),0.1)
#plot(range(1,Nexpansion+1), AveragePhase[discard:len(AveragePhase)])
#plt.hold('True')
#xlabel('$N$')
#ylabel('$\langle \phi \\rangle$')
#plot(xinterp, [Exponential(xinterp[i] ,aPhase, bPhase, cPhase) for 
#               i in range(len(xinterp))])
#plt.figure(2)
#plot(range(1,Nexpansion+1), AveragePeriod[discard:len(AveragePeriod)])
#plt.hold('True')
#xlabel('$N$')
#ylabel('$\langle T \\rangle$')
#plot(xinterp, [Exponential(xinterp[i] ,aPeriod, bPeriod, cPeriod) for 
#               i in range(len(xinterp))])
plt.show()
