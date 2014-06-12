"""
    Main file to produce results that are going to be included in the 2modes paper
    Not every calculation needs to be done at each run, data can be produced
    and stored locally, see the booleans in line 18
    
    Requires numpy ver > 1.8
"""

import numpy as np
from scipy import interpolate, integrate
from scipy.misc import derivative
from scipy.optimize import newton, fsolve, fmin
import sys
import sqlite3
#Initiate plotting environment:
import matplotlib as mpl
from pylab import plot, xlabel, ylabel, show, savefig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from subprocess import call
#Import twomode system module
import twomode 
from PADS import Lyndon
import oct2py
import os
peddir = os.getcwd()+'/ped' #Add /ped to the current working directory location
oct2py.octave.addpath(peddir) #Add mfiles for ped

#Booleans:
computeSolution = False
computePsect = False
computeArcLengths = False
computeRPO = True
computeRPOred3 = False
plotPsect = False
plotRetmap = False

#Search parameters:
nPrimeMax = 8 #Will search for [1,m]-cycles

#Only relative equilibrium:
reqv = np.array([0.43996557973671596,
                 0,
                 -0.38626705952942764,
                 0.07020440068369532], float) 
reqv3D = twomode.four2three(reqv)
np.savetxt('data/reqv3d.dat', reqv3D)
#Compute reduced stability matrix for the relative equilibrium:
Ared = twomode.StabilityMatrixRed(reqv)
#Compute eigenvalues and eigenvectors:
w, v = np.linalg.eig(Ared)
#Pick the real part of unstable eigenvector as the Poincare section direction:
unstabledir2 = np.real(v[:,0])
unstabledir = np.imag(v[:,0])
#unstabledir = np.real(w[0]*v[:,0]+w[1]*v[:,1])
unstabledir = unstabledir/np.linalg.norm(unstabledir) #Normalize
unstabledir2 = unstabledir2/np.linalg.norm(unstabledir2)
#Auxiliary vector to define Poincare basis:
vaux3D = np.array([0,0,1], float)
unstabledir3D = twomode.four2three(unstabledir)
unstabledir23D = twomode.four2three(unstabledir2)
np.savetxt('data/unstabledir3d.dat', unstabledir3D)
np.savetxt('data/unstabledir23d.dat', unstabledir23D)
nhat3D = np.cross(vaux3D, unstabledir3D)
nhat = twomode.three2four(nhat3D)
vaux = twomode.three2four(vaux3D)
#Poincare Basis:
psbasis = np.array([unstabledir, vaux, nhat]).transpose()
#Lie algebra generator:
T = twomode.generator()
xp = np.array([1,0,0,0], float)
tp = np.dot(T,xp)

if computeSolution:
    
    tf = 2000;
    dt = 0.01;
    epsilon = 1e-2;
    x0 = reqv+epsilon*unstabledir
    #x0 = np.array([1e-3,0,1e-3,0], float)
    xphi0 = np.append(x0, 0)
    print xphi0
    
    t = np.linspace(0, tf, np.floor(tf/dt)+1)
    xphisol = twomode.intslice(xphi0, t)

    txphisol = np.concatenate((t.reshape(-1,1),xphisol), axis=1)
    txphisol = txphisol[20000:-1,:]
    #Create a figure window
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #Modify axis colors:
    ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
    ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
    ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
            
    ax.plot(xphisol[:,0], 
    xphisol[:,2], 
    xphisol[:,3], linewidth=0.2, color='#3c5f96')
    
    ax.set_xlabel('\n $\hat{x}_1$ \t  ', fontsize=32)
    ax.set_ylabel('\n $\hat{x}_2$ \t', fontsize=32)
    ax.set_zlabel('$\hat{y}_2$   ', fontsize=32)
    
    np.savetxt('data/txphisol.dat', txphisol)
    
    plt.show()
    
if computePsect:
    import poincare
    if not('txphisol' in locals()):
        txphisol = np.loadtxt('data/txphisol.dat')
    ps = poincare.computeps(txphisol, reqv, nhat, 1)
    np.savetxt('data/ps.dat', ps)

if not('ps' in locals()):
    ps = np.loadtxt('data/ps.dat')

#Project onto poincare basis:
ps2D = np.dot(psbasis.transpose(), (ps[:,1:5]-reqv).transpose()).transpose()
ps2D = ps2D[:,0:2]
#Sort ps2D for interpolation:
sortingindices = np.argsort(ps2D[:,0])
ps2Dsorted = np.array(ps2D[sortingindices, :],float)

#Interpolate to psect:
print 'Interpolating the Poincare section'
tckps = interpolate.splrep(ps2Dsorted[:,0],ps2Dsorted[:,1], k=3)
dxintps = (ps2Dsorted[np.size(ps2Dsorted,0)-1,0] - ps2Dsorted[0,0])/10000
xintps = np.arange(ps2Dsorted[0,0], 
                      ps2Dsorted[np.size(ps2Dsorted,0)-1,0]+dxintps, 
                      dxintps)
yintps = interpolate.splev(xintps, tckps)

def psarclength(x):
    """
    Function to compute arclenghts in the Poincare section
    Ref: http://www.mathwords.com/a/arc_length_of_a_curve.htm 
    """
    def ds(x):
        dypsect = derivative(interpolate.splev, x, dx=1e-6, args=(tckps,), order=5)
        return np.sqrt(1 + dypsect**2)
    #quad, err = integrate.quad(ds, ps2Dsorted[0,0], x)
    return x
    #return quad

if computeArcLengths:
    print 'Computing arclengths corresponding to data'
    #Find the first data point for the Arclengths to discard the transients
    iArcLength0 = np.argwhere(ps[:,0]>200)[0]
    sn = np.array([psarclength(ps2D[i,0]) for i in 
    range(iArcLength0, np.size(ps2D,0))], float)
    snmin = np.min(sn)
    snmax = np.max(sn) - snmin
    #sn = np.array([sn[i] for i in range(np.size(sn, 0))])
    snplus1 = sn[1:]
    sn = sn[0:-1]
    RetMapData = np.array([sn, snplus1], float).transpose()
    np.savetxt('data/RetMapData.dat', RetMapData)
    snMinMax = [snmin, snmax]
    np.savetxt('data/RetMapMinMax.dat', snMinMax)

if not('sn' in locals()):
    RetMapData = np.loadtxt('data/RetMapData.dat')
    snmin, snmax = np.loadtxt('data/RetMapMinMax.dat')
    sn = RetMapData[:,0]
    snplus1 = RetMapData[:,1]

print 'Interpolating to the return map'
isortRetMap = np.argsort(sn)
tckRetMap = interpolate.splrep(sn[isortRetMap],snplus1[isortRetMap], k=1)
dxintRetMap = (np.max(sn)-np.min(sn))/100000
xintRetMap = np.arange(np.min(sn), np.max(sn), dxintRetMap)
yintRetMap = interpolate.splev(xintRetMap, tckRetMap)

#Return map function:
def retmap(sn):
    snp1 = interpolate.splev(sn, tckRetMap)
    return float(snp1)
#nth return map function:
def retmapm(n, sn):
    if n == 0:
        return sn
    else:
        snpn = retmap(sn)
        for i in range(n - 1):
            snpn = retmap(snpn)
        return  snpn

print "Computing the kneading sequence"
nMax = nPrimeMax;
sCritical = fmin(lambda x: -interpolate.splev(x, tckRetMap), 1)
print "sCritical:"
print sCritical
def fCritical(s):
    po = retmapm(3, s) - s
    return po
s3Critical = newton(fCritical, sCritical*0.999, tol=1.48e-12)
#s3Critical = newton(fCritical, sCritical*1.001, tol=1.48e-12)
#s3Critical = sCritical*0.9995
print "s3Critical:"
print s3Critical

Kneading = np.copy(s3Critical)
KneadingSequence = []
KneadingValueBin = '0.'
KneadingValue = 0
for i in range(nMax):
    xnext = retmap(Kneading[-1])
    Kneading=np.append(Kneading, xnext)
    if xnext > sCritical:
        KneadingSequence.append(1)
        if i == 0:
            KneadingValueBin = KneadingValueBin+'1'
        else:
            KneadingValueBin = KneadingValueBin+str(int(not(int(KneadingValueBin[-1]))))
    else:
        KneadingSequence.append(0)
        if i == 0:
            KneadingValueBin = KneadingValueBin+'0'
        else:
            KneadingValueBin = KneadingValueBin+KneadingValueBin[-1]
            
    KneadingValue = KneadingValue + (int(KneadingValueBin[-1])*0.5)**(i+1)
print "Kneading:"
print Kneading
print "Kneading Sequence:"
print KneadingSequence
print "Kneading Value (binary):"
print KneadingValueBin
print "Kneading Value (real):"
print KneadingValue

def Itinerary(s, n):
    """
    Compute future itinerary of a point s on the return map
    """
    itinerary = []
    for i in range(n):
        s = retmap(s)
        if s>sCritical:
            itinerary.append(1)
        elif s<sCritical:
            itinerary.append(0)
    return itinerary

def Splus2gamma(itinerary):
    gamma = 0
    for i in range(len(itinerary)):
        if i == 0:
            gammaBin = '0.'+str(itinerary[i])
        elif itinerary[i]==0:
            gammaBin = gammaBin + gammaBin[-1]
        elif itinerary[i]==1:
            gammaBin = gammaBin + str(int(not(int(gammaBin[-1]))))      
        gamma = gamma + (int(gammaBin[-1])*0.5)**(i+1)
    return gamma, gammaBin
    
def TopologicalCoordinate(x, n):
    """
    Compute topological coordinate of a point in the return map
    """
    itinerary = Itinerary(x,n)
    gamma, gammaBin = Splus2gamma(itinerary)
    return gamma

KneadingValue = TopologicalCoordinate(s3Critical, nPrimeMax)
print "Kneading value fun:"
print KneadingValue

PrimeCycles = []

for j in range(1,nPrimeMax+1):
    for lyndon in Lyndon.LyndonWordsWithLength(2,j):
        PrimeCycles.append([j, list(lyndon) ])

nRPO = -1
AdmissibleCycles = []
#Compute the maximum topological coordinate for each prime cycle:
for k in range(len(PrimeCycles)):
    MaxTopCoord = 0
    nCycle = len(PrimeCycles[k][1])
    for j in range(nCycle):
        Permutation = [PrimeCycles[k][1][i - j] for i in range(nCycle)]
        TopCoorPerm = Splus2gamma((Permutation*(nPrimeMax/nCycle+1))[0:nPrimeMax])[0]
        if TopCoorPerm > MaxTopCoord:
            MaxTopCoord = TopCoorPerm
    PrimeCycles[k].append(MaxTopCoord)
    if MaxTopCoord > KneadingValue:
        PrimeCycles[k].append('Inadmissible')
        nRPO = nRPO + 1
    else:
        PrimeCycles[k].append('Admissible')
        AdmissibleCycles.append([PrimeCycles[k][0], PrimeCycles[k][1]])
AdmissibleCycles = AdmissibleCycles[1:]

orig_stdout = sys.stdout

f = file('data/primeCycles.txt', 'w')
sys.stdout = f

print "PrimeCycles upto length "+str(nPrimeMax)
for i in range(len(PrimeCycles)):
    print PrimeCycles[i]

sys.stdout = orig_stdout
f.close()

def fs2ps2D(px, s):
    """
        Function to be solved to get relative x coordinate of the rpo on
        the Poincare section
    """
    sfun = psarclength(px)  
    return sfun-s

#From arclength to projected Poincare section coordinates:
def s2ps2D(s):
    #Guessing the initial point:
    i0 = np.argmin(np.absolute(sn - s))
    #p0x = ps2D[i0,0]
    p0x = s

    px = newton(fs2ps2D, p0x, args=(s,), tol=1e-12)
    py = interpolate.splev(px, tckps)
    
    return np.array([px, py])

def ps2D2psxhat(ps2Di):
        """
            Takes the relative position projected on Poincare section basis 
            on the Poincare section and returns the position on the reduced 
            state space
        """
        psrelproj3D = np.array([ps2Di[0], ps2Di[1], 0], float)
        psrel = np.dot(psbasis, psrelproj3D)
        ps = psrel + reqv
        #ps = np.array([ps3D[0], 0, ps3D[1], ps3D[2]], float)
        return ps

def s2xhat(s):
    """
    From arclength to reduced state space coordinates.
    """
    ps2D = s2ps2D(s)
    xhat = ps2D2psxhat(ps2D)    
    return xhat

def timeofflight(xhat0):
    """
    Computes time of flight for xhat0 on the Poincare section to come back onto
    the section for a second time 
    """
    import poincare
    #Find the point on the computed poincare section, closest to the x0
    imin = np.argmin(np.linalg.norm(ps[:,1:5]-xhat0, axis=1))
    #Take its time of flight as
    if imin < np.size(ps,0)-1: 
        Tapproximate = ps[imin+1,0]-ps[imin,0]
    else:
        Tapproximate = ps[imin,0]-ps[imin-1,0]
    #Integrate for a little bit longer than the approximated integration time:
    stoptime = 1.2*Tapproximate
    numpoints = int(stoptime/0.01)
    #Integration time array:
    t = np.linspace(0, stoptime, numpoints)
    xhatphi0 = np.append(xhat0, np.array([0], float))
    xhatsol = twomode.intslice(xhatphi0, t, abserror=1.0e-14, relerror=1.0e-12)
    tx = np.append(np.array([t], float).transpose(), xhatsol, axis=1)
    #Compute Poincare section:
    psreturn=poincare.computeps(tx, reqv, nhat, 1)
    #print psreturn
    #Take the time nearest to the approximated time. This is due to the
    #fact that the array ps sometimes includes the initial point and sometimes
    #does not, hence we are not always sure the position of the first return.
    itof = np.argmin(np.abs(psreturn[:,0]-Tapproximate))
    tof = psreturn[itof,0]
    return tof
    
def phireturn(xhat0, tof):
        """
        Computes phi for xhat0 on the Poincare section to come back onto
        the section for a second time 
        """
        stoptime = tof
        numpoints = 2
        #Integration time array:
        t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]        
        xsol = twomode.intfull(xhat0, t, abserror=1.0e-14, relerror=1.0e-12)
        #Phase of the first mode is the slice phase
        phi = np.angle(xsol[1,0] + 1j*xsol[1,1])    
        return -phi        

Adaptive = True
factor = 2

if computeRPO:
    smin = np.min(sn)
    smax = np.max(sn)
    #scandidates = np.zeros([1,2]) #dummy matrix to hold po candidate arclengths
    scandidates = [] #dummy matrix to hold po candidate arclengths
    
    for i in range(nPrimeMax):
        #Define the function zeros of which would correspond to the periodic 
        #orbit arclengths:
        def fpo(s):
            po = retmapm(i+1, s) - s
            return po
            
        fpoevo=0
        for s0 in np.arange(smin, smax, (smax-smin)/40000):
            fpoev = fpo(s0)
            if fpoev * fpoevo < 0: #If there is a zero-crossing, look for the root: 
                sc = newton(fpo, s0, tol=1.48e-8)
                #print "sc = %f" %sc
                newcandidate = 1
                for j in range(len(scandidates)):
                    #Discard if found candidate is previously found:
                    if np.abs(scandidates[j][1] - sc)<1e-9: 
                        newcandidate=0
                if newcandidate:
                    CandidateItinerary = Itinerary(sc, i+1)
                    scandidates.append([i+1, sc, Itinerary(sc, i+1)])
                    for k in range(len(AdmissibleCycles)):
                        if CandidateItinerary == AdmissibleCycles[k][1]:
                            AdmissibleCycles[k].append([retmapm(n, sc) for n 
                                                       in range(i+1)])
                            AdmissibleCycles[k].append(np.array([s2xhat(
                             AdmissibleCycles[k][2][l]) for l in range(i+1)]))
                            AdmissibleCycles[k].append(np.array([timeofflight(
                             AdmissibleCycles[k][3][l]) for l in range(i+1)],
                                                                        float))
                            AdmissibleCycles[k].append(np.array([phireturn(
                             AdmissibleCycles[k][3][l],
                             AdmissibleCycles[k][4][l]) for l in range(i+1)],
                                                                        float))
                    #scandidates = np.append(scandidates, np.array([[i+1, sc]], float), axis=0)
                    
            fpoevo = fpoev        


    print "Admissible Cycles upto length before subdivisions:"
    for i in range(len(AdmissibleCycles)):
        print AdmissibleCycles[i]
        # #Divide intervals into smaller subintervals for multiple shooting:
    #nsub = 1 #number of subintervals
    #for k in range(len(AdmissibleCycles)):
        #l = 0
        #while l < len(AdmissibleCycles[k][3]):
            ##tau:
            #AdmissibleCycles[k][4][l] = AdmissibleCycles[k][4][l]/float(nsub)
            ##phi:
            #AdmissibleCycles[k][5][l] = phireturn(AdmissibleCycles[k][3][l],
                                                  #AdmissibleCycles[k][4][l])
            #for m in range(nsub-1):
                ##x:
                #AdmissibleCycles[k][3] = \
                #np.insert(AdmissibleCycles[k][3], l+m+1,
                          #twomode.ftauRed(AdmissibleCycles[k][3][l+m,:],
                                #AdmissibleCycles[k][4][l]),
                          #axis=0)
                ##tau:
                #AdmissibleCycles[k][4] = np.insert(AdmissibleCycles[k][4],
                                                   #l+m+1,
                                                   #AdmissibleCycles[k][4][l])
                ##phi:
                #AdmissibleCycles[k][5] = np.insert(AdmissibleCycles[k][5],
                                    #l+m+1,
                                    #phireturn(AdmissibleCycles[k][3][l+m+1],
                                    #AdmissibleCycles[k][4][l+m+1]))
            #l += nsub

    print "Admissible Cycles upto length "+str(nPrimeMax)
    for i in range(len(AdmissibleCycles)):
        print AdmissibleCycles[i]

    #raw_input('Press enter to continue ...')

    print "Starting the Newton search..."
    tol = 1e-9
    #for i in range(1,2):
    for i in range(len(AdmissibleCycles)):

        converged = False

        x = AdmissibleCycles[i][3]
        tau = AdmissibleCycles[i][4]
        phi = AdmissibleCycles[i][5]
        nCycle = len(x)
        #Error vector
        while not(converged):
            Error = np.array([np.append(x[(k+1)%nCycle] - np.dot(twomode.LieElement(phi[(k)%nCycle]),
                                                               twomode.ftau(x[(k)%nCycle], tau[(k)%nCycle])),
                             np.array([0, 0], float)) for k in range(nCycle)], float)

            Error = Error.reshape(np.size(Error))
            print "Error"
            print Error
            #raw_input("Press enter to continue...")
            if np.max(np.abs(Error)) < tol:
                converged = True

            N = np.size(Error,0)
            nDim = np.size(x[0], 0)
            #A-Matrix:
            A = np.zeros((N,N))
            for k in range(nCycle):
                A[(nDim+2)*k: (nDim+2)*k + nDim, (nDim+2)*k: (nDim+2)*k + nDim] = np.dot(twomode.LieElement(phi[k]),
                                                                                        twomode.Jacobian(x[k], tau[k]))
                A[(nDim+2)*k: (nDim+2)*k + nDim, (nDim+2)*k + nDim] = np.dot(twomode.LieElement(phi[k]),
                                                                        twomode.vfullssp(twomode.ftau(x[k], tau[k])))
                A[(nDim+2)*k: (nDim+2)*k + nDim, (nDim+2)*k + nDim + 1] = np.dot(twomode.generator(),
                                                                                 np.dot(twomode.LieElement(phi[k]),
                                                                                 twomode.ftau(x[k], tau[k])))
                A[(nDim+2)*k + nDim, (nDim+2)*k: (nDim+2)*k + nDim] = twomode.vfullssp(x[k])
                A[(nDim+2)*k + nDim + 1, (nDim+2)*k: (nDim+2)*k + nDim] = np.dot(twomode.generator(), x[k])
                #A[(nDim+2)*k + nDim, (nDim+2)*k: (nDim+2)*k + nDim] = nhat
                #A[(nDim+2)*k + nDim + 1, (nDim+2)*k: (nDim+2)*k + nDim] = np.array(tp)

                A[(nDim+2)*((k)%nCycle): (nDim+2)*((k)%nCycle) + nDim,
                  (nDim+2)*((k+1)%nCycle): (nDim+2)*((k+1)%nCycle) + nDim] = \
                A[(nDim+2)*((k)%nCycle): (nDim+2)*((k)%nCycle) + nDim,
                  (nDim+2)*((k+1)%nCycle): (nDim+2)*((k+1)%nCycle) + nDim] - np.identity(nDim)
    
            #Compute Deltas:
            Delta=np.dot(np.linalg.inv(A), Error)
            print "Delta"
            print Delta
            
            converging = False
            if Adaptive:
                xx = np.empty(np.shape(x))
                tautau = np.zeros(np.shape(tau))
                phiphi = np.zeros(np.shape(phi))
                iAdaptive = 0
                while not(converging):
                    iAdaptive = iAdaptive + 1
                    print "iAdaptive:"
                    print iAdaptive
                    for k in range(nCycle):
                        xx[k] = x[k] + Delta[(nDim+2)*k:(nDim+2)*k+nDim]
                        tautau[k] = tau[k] + Delta[(nDim+2)*k+nDim]
                        phiphi[k] = phi[k] + Delta[(nDim+2)*k+nDim+1]
                    
                    ErrorNext = np.array([xx[(k+1)%nCycle] - \
                    np.dot(twomode.LieElement(phiphi[(k)%nCycle]), \
                    twomode.ftau(xx[(k)%nCycle], tautau[(k)%nCycle])) \
                    for k in range(nCycle)], float)

                    
                    if np.max(np.abs(ErrorNext)) < np.max(np.abs(Error)):
                        converging = True
                    else:
                        Delta = Delta * (factor**(-1))
                        

            #Update:
            for k in range(nCycle):
                x[k] = x[k] + Delta[(nDim+2)*k:(nDim+2)*k+nDim]
                tau[k] = tau[k] + Delta[(nDim+2)*k+nDim]
                phi[k] = phi[k] + Delta[(nDim+2)*k+nDim+1]

    #np.savetxt('data/AdmissibleCycles.dat', AdmissibleCycles)

    #Create a database and write RPOs in it:
    
    conn = sqlite3.connect('data/rpo.db')
    c = conn.cursor()
    c.execute("DROP TABLE IF EXISTS rpos")
    c.execute(" CREATE TABLE rpos (rpono int, itinerary text, x1 real, \
                y1 real, x2 real, y2 real, period real, phase real, \
                FloquetMults text) ")
                #FloquetMults text, FloquetExps text) ")

    for i in range(len(AdmissibleCycles)):
        
        rpo = AdmissibleCycles[i][3][0]
        Period = sum(AdmissibleCycles[i][4])
        Phase = sum(AdmissibleCycles[i][5])
        JJ = np.dot(twomode.LieElement(AdmissibleCycles[i][5][0]),
                    twomode.Jacobian(AdmissibleCycles[i][3][0], 
                                     AdmissibleCycles[i][4][0]))
        for k in range(1, len(AdmissibleCycles[i][3])):
            JJ = np.append(JJ, 
                          np.dot(twomode.LieElement(AdmissibleCycles[i][5][k]),
                          twomode.Jacobian(AdmissibleCycles[i][3][k], 
                                           AdmissibleCycles[i][4][k])), axis=1)
        #raw_input("dfasfdsavfass")
        eig = oct2py.octave.ped(JJ)
        FloquetMults = np.array([np.real(np.exp(eig[j, 0])*eig[j, 1]) for j in range(len(eig))], 
                                float)
        #FloquetMults = twomode.Floquet(rpo, Period, Phase)
        print FloquetMults
        #FloquetExps = np.log(FloquetMults)

        itinerary = ''.join(map(str, AdmissibleCycles[i][1]))
        c.execute("INSERT INTO rpos VALUES ("+str(i+1)+", '" \
                   +itinerary+"', " \
                   +', '.join(map(str, rpo))+", " \
                   +str(Period)+", " \
                   +str(Phase)+", '" \
                   +str(FloquetMults)+"')")
                   #+str(FloquetMults)+"' , '" \
                   #+str(FloquetExps)+"' )")

        c.execute("DROP TABLE IF EXISTS rpo"+itinerary)
        c.execute(" CREATE TABLE rpo"+itinerary+" (x1 real, y1 real, x2 real, \
        y2 real, tau real, phi real)")
    
        for k in range(len(AdmissibleCycles[i][3])):
            query = "INSERT INTO rpo"+itinerary+" VALUES("+  \
                      ', '.join(map(str, AdmissibleCycles[i][3][k]))+", " \
                      +str(AdmissibleCycles[i][4][k])+", " \
                      +str(AdmissibleCycles[i][5][k])+")"              
            c.execute(query)
    
    conn.commit()
    conn.close()

if plotPsect:
    fig=plt.figure(1, figsize=(8,8))
    plot(ps2D[:,0], ps2D[:,1], '.b', ms=10)
    plt.hold(True)
    plot(xintps, yintps, 'k', lw=2)
    ax=fig.gca()
    #ax.set_aspect('equal')
    ax.set_xlabel('$v_1$', fontsize=24)
    ax.set_ylabel('$v_{2 \quad (\\times 100)}$', fontsize=24)
    Nticks=5
    #xticks = np.linspace(min(ps2D[:,0]), max(ps2D[:,0]), Nticks)
    #ax.set_xticks(xticks)
    #ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=16);
    #yticks = np.linspace(min(ps2D[:,1]), max(ps2D[:,1]), Nticks)
    #ax.set_yticks(yticks)
    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #yticks=yticks*100
    #ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=16);

    savefig('Psect.pdf', bbox_inches='tight', dpi=100)
    call(["pdfcrop", "Psect.pdf", "Psect.pdf"], shell=True)

    show()

if plotRetmap:
    fig=plt.figure(1, figsize=(8,8))
    srange = np.arange(np.min(sn), np.max(sn), np.max(sn)/50000)
    plot(sn, snplus1, '.b',                                                     ms=10)
    plt.hold(True)
    plot(xintRetMap, yintRetMap, 'k', lw=2)
    plot(srange, srange, 'g', lw=2)

    nKneading = np.size(Kneading)
    for i in range(nKneading-1):
        pair1 = [Kneading[i], Kneading[(i+1)]]
        pair2 = [Kneading[i], Kneading[(i+1)]]
        plot(np.linspace(min(pair1), max(pair1), 10),
            [pair1[1] for k in range(10)], '--r', lw=1.5)
        plot([pair2[0] for k in range(10)],
            np.linspace(min(pair2), max(pair2), 10), '--r', lw=1.5)
    
    #plot(srange, [0.825 for  i in range(np.size(srange,0))], 'r')
    ax = fig.gca()
    #ax.set_aspect('equal')
    smin = np.min(sn)
    smax = np.max(sn)
    ax.set_xlim(smin,smax)
    ax.set_ylim(smin,smax)
    ax.set_xlabel('$s_n$', fontsize=24)
    ax.set_ylabel('$s_{n+1}$', fontsize=24)
    #Nticks = 5

    #xticks = np.linspace(smin, smax, Nticks)
    #ax.set_xticks(xticks)
    #ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=16); 

    #yticks = np.linspace(smin, smax, Nticks)
    #ax.set_yticks(yticks)
    #ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=16); 
    
    #plt.figure(2, figsize=(8,8))
    #sp3 = np.array([retmapm(7, sn) for sn in srange])
    #plot(srange, sp3, 'b')
    #plt.hold(True)
    #plot(srange,srange,'g')
    
    #savefig('RetMap.pdf', bbox_inches='tight', dpi=100) 
    #call(["pdfcrop", "RetMap.pdf", "RetMap.pdf"], shell=True)
    
    show()
