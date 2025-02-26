import numpy as np
import sqlite3
from StringIO import StringIO 
from subprocess import call
from scipy.optimize import fsolve
import sys
import sympy
#Initiate plotting environment:
import matplotlib as mpl
from pylab import plot, xlabel, ylabel, show, savefig
import matplotlib.pyplot as plt
import twomode

FiniteGrammar=False

if FiniteGrammar:
    Ncycle = 79
    NmaxExp = 12
    conn = sqlite3.connect('data/rpo.db')

else:
    Ncycle = 84
    NmaxExp = 11
    conn = sqlite3.connect('data/rpoall.db')

c = conn.cursor()
rpos = []

for rpono in range(1,Ncycle + 1):
    c.execute("SELECT * FROM rpos WHERE rpono = "+str(rpono))
    a = c.fetchall()
    
    rpono,  itinerary, x1, y1, x2, y2, T, phifetch , floquetFetch= a[0]
    #Format the Floquet exponents such that it would be readable by loadtxt as
    #a numpy array
    floquetobj = StringIO(str(floquetFetch).strip('[ ]'))
    floquet = np.loadtxt(floquetobj)
    TopLength = len(str(itinerary))
    
    c.execute("SELECT phi FROM rpo"+str(itinerary))
    a = c.fetchall()
    phi = -sum([float(str(a[i]).strip('( ,)')) for i in range(len(a))])
    print rpono, phi, T
    rpos.append([TopLength, T, floquet[0], phi])

conn.close()

s, z, beta = sympy.symbols('s z beta')
EscapeRate = []
AveragePhase = []
AveragePhaseVar = []
AveragePeriod = []
AveragePhaseSpeed = []
AverageLyapunovExponent = []
ConservationRule = []

for Nexpansion in range(1,NmaxExp+1):
#for Nexpansion in range(NmaxExp,NmaxExp+1):
    SpectralDeterminant = 1
    AvgT = 1
    Avgphi = 1
    AvgLy = 1
    #Nexpansion = 3
    Exponent = 0
    #Compute Spectral Determinant according to ChaosBook v14 eq. 19.6
    for i in range(Ncycle):
        npr = rpos[i][0]
        Tp = rpos[i][1]
        Lambda = rpos[i][2]
        phi = rpos[i][3]
        Sum = 0
        SumT = 0
        Sumphi = 0
        SumLy = 0
        r = 1
        #if npr == 3:
        #    continue
        
        while npr*r <= Nexpansion:
            Sum = Sum + (sympy.exp(-r*s*Tp)*z**(npr*r))/(r*abs(1.0-Lambda**r))
            #SumT = Sum + (sympy.exp(r*(beta*Tp-s*Tp))*z**(npr*r))/(r*abs(1.0-Lambda**r))
            Sumphi = Sumphi + (sympy.exp(r*(beta*phi - s*Tp))*z**(npr*r))/(r*abs(1.0-Lambda**r))
            SumLy = SumLy + (sympy.exp(r*(beta*(np.log(np.abs(Lambda))) - s*Tp))*z**(npr*r))/(r*abs(1.0-Lambda**r))
            r += 1
        #Expand the exponentials and discard higher order terms:
        print "Exponentiating sums"
        ExpSum = sympy.series(sympy.exp(-Sum), z, n=Nexpansion+1).subs(sympy.O(z**(Nexpansion+1)), 0)
        print "done with first"
        #ExpSumT = sympy.series(sympy.exp(-SumT), z, n=Nexpansion+1).subs(sympy.O(z**(Nexpansion+1)), 0)
        ExpSumphi = sympy.series(sympy.exp(-Sumphi), z, n=Nexpansion+1).subs(sympy.O(z**(Nexpansion+1)), 0)
        print "done with second"
        ExpSumLy = sympy.series(sympy.exp(-SumLy), z, n=Nexpansion+1).subs(sympy.O(z**(Nexpansion+1)), 0)
        print "Expanding exponentials"
        SpectralDeterminant = (SpectralDeterminant * ExpSum).expand()
        #AvgT = (AvgT * ExpSumT).expand()
        Avgphi = (Avgphi * ExpSumphi).expand()
        AvgLy = (AvgLy * ExpSumLy).expand()
        print "Fixing degrees"
        while sympy.degree(SpectralDeterminant, z) > Nexpansion:
            SpectralDeterminant = SpectralDeterminant - sympy.LT(SpectralDeterminant, z)
            SpectralDeterminant = sympy.collect(SpectralDeterminant, z)
    
        #while sympy.degree(AvgT, z) > Nexpansion:
        #    AvgT = AvgT - sympy.LT(AvgT, z)
        #    AvgT = sympy.collect(AvgT, z)

        while sympy.degree(Avgphi, z) > Nexpansion:
            Avgphi = Avgphi - sympy.LT(Avgphi, z)
            Avgphi = sympy.collect(Avgphi, z)
    
        while sympy.degree(AvgLy, z) > Nexpansion:
            AvgLy = AvgLy - sympy.LT(AvgLy, z)
            AvgLy = sympy.collect(AvgLy, z)
    
    SpectralDeterminant = SpectralDeterminant.subs(z,1)    
    #AvgT = AvgT.subs(z,1)    
    Avgphi = Avgphi.subs(z,1)    
    AvgLy = AvgLy.subs(z,1)    
        #Exponent = Exponent - Sum
    #Extract coefficients of the expansion:
    #SpectralDeterminant = sympy.series(sympy.exp(Exponent),z,n=Nexpansion+1).subs(sympy.O(z**(Nexpansion+1)), 0)
    #SpectralDeterminant = sympy.series(sympy.exp(Exponent),z,n=Nexpansion+1).subs(sympy.O(z**(Nexpansion+1)), 0).subs(z, 1)
    ConservationRule.append(SpectralDeterminant.subs(s,0))
    #EscapeRate.append(float(sympy.nsolve(SpectralDeterminant, 0)))
    f = sympy.lambdify(s, SpectralDeterminant, "numpy")
    def fcomplex(xy):
        x, y = xy
        fcomp = np.array([np.real(f(x + 1j*y)), np.imag(f(x + 1j*y))])
        return fcomp
    #found = False
    #step = 1e-4
    #splusnext = 0
    #sminusnext = 0
    #while not(found):
        #splus = splusnext
        #sminus = sminusnext
        #splusnext += step
        #sminusnext -= step
        #if f(sminus)*f(sminusnext) < 0:
            #EscapeRate.append(float(fsolve(f, (sminus + sminusnext)/2.0)))
            #found = True
        #if f(splus)*f(splusnext) < 0:
            #EscapeRate.append(float(fsolve(f, (splus + splusnext)/2.0)))
            #found = True
    EscapeRate.append(float(fsolve(f, 0)))
    #AvgT = sympy.diff(SpectralDeterminant, s).subs(s, EscapeRate[-1])
    AveragePeriod.append(sympy.diff(SpectralDeterminant, s).subs(s, 
                        EscapeRate[-1]))
    #AveragePhase.append(Avgphi.subs(s, EscapeRate[-1]))
    AveragePhase.append(sympy.diff(-Avgphi, beta).subs(beta,0).subs(s,
                        EscapeRate[-1]))
    AveragePhaseVar.append(sympy.diff(-Avgphi, beta, 2).subs(beta,0).subs(s, 
                           EscapeRate[-1]))
    AverageLyapunovExponent.append((sympy.diff(-AvgLy, beta).subs(beta,0).subs(s,
                           EscapeRate[-1]))/AveragePeriod[-1])
    #raw_input("stop!")
    AveragePhaseSpeed.append(AveragePhase[-1]/AveragePeriod[-1])
    #EscapeRate.append(fsolve(fcomplex, [0, 0]))
    print "EscapeRate", EscapeRate
    print "AveragePeriod", AveragePeriod
    print "AveragePhase", AveragePhase
    print "AveragePhaseVar", AveragePhaseVar
    print "AveragePhaseSpeed", AveragePhaseSpeed
    print "AverageLyapunovExponent", AverageLyapunovExponent

#f = open("tex/s0.tex", "w")

#f.write("\\begin{table}\n")
#f.write("\t\\begin{tabular}{c|c}\n")

#f.write("\t N & s_0 \\\\ \n")
#f.write("\t\\hline\n")

#for i in range(len(EscapeRate)):
    #f.write("\t%s & " % str(i+1))
    #f.write("%5.8f \\\\ \n " % float(EscapeRate[i]))
    
#f.write("\t\\end{tabular}\n")
#f.write("\t\\caption{Leading zero of the spectral determinant (\\beta = 0) \
#computed using the finite grammar approximation}\n")
#f.write("\t\\label{t-s0}\n")
#f.write("\\end{table}")

#f.close()     

    DiffusionCoefficient = [AveragePhaseVar[i]/(2*AveragePeriod[i]) for i in range(len(AveragePhaseVar))]
    
    escrate = np.array(EscapeRate, float)
    escrate = -escrate
    np.savetxt('data/EscapeRate.dat', escrate)
    np.savetxt('data/AveragePeriod.dat', AveragePeriod)
    np.savetxt('data/AveragePhase.dat', AveragePhase)
    np.savetxt('data/AveragePhaseVar.dat', AveragePhaseVar)
    np.savetxt('data/AverageLyapunovExponent.dat', AverageLyapunovExponent)
    np.savetxt('data/DiffusionCoefficient.dat', DiffusionCoefficient)
    np.savetxt('data/AveragePhaseSpeed.dat', AveragePhaseSpeed)
    
