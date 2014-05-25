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
#Initiate plotting environment:
import matplotlib as mpl
from pylab import plot, xlabel, ylabel, show, savefig
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from subprocess import call
#Import twomode system module
import twomode 

#Booleans:
computeSolution = False
computePsect = False
computeArcLengths = False
computeRPO = False
plotPsect = False
plotRetmap = False

#Search parameters:
m = 1 #Will search for [1,m]-cycles

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
unstabledir = np.real(v[:,0])
unstabledir2 = np.imag(v[:,0])
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
    
    tf = 1000;
    dt = 0.01;
    epsilon = 1e-3;
    x0 = reqv+epsilon*unstabledir
    xphi0 = np.append(x0, 0)
    print xphi0
    
    t = np.linspace(0, tf, np.floor(tf/dt)+1)
    xphisol = twomode.intslice(xphi0, t)

    txphisol = np.concatenate((t.reshape(-1,1),xphisol), axis=1)
    #Create a figure window
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #Modify axis colors:
    ax.w_xaxis.set_pane_color((1, 1, 1, 1.0))
    ax.w_yaxis.set_pane_color((1, 1, 1, 1.0))
    ax.w_zaxis.set_pane_color((1, 1, 1, 1.0))
            
    ax.plot(xphisol[:,0], 
    xphisol[:,2], 
    xphisol[:,3], linewidth=0.1, color='#3c5f96')
    
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
tckps = interpolate.splrep(ps2Dsorted[:,0],ps2Dsorted[:,1])
dxintps = (ps2Dsorted[np.size(ps2Dsorted,0)-1,0] - ps2Dsorted[0,0])/100
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
    quad, err = integrate.quad(ds, ps2Dsorted[0,0], x)  
    return quad

if computeArcLengths:
    print 'Computing arclengths corresponding to data'
    #Find the first data point for the Arclengths to discard the transients
    iArcLength0 = np.argwhere(ps[:,0]>300)[0]
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
    return snp1
#nth return map function:
def retmapm(n, sn):
    snpn = retmap(sn)
    for i in range(n - 1):
        snpn = retmap(snpn)
    return  snpn

print "Computing the kneading sequence"
nMax = 3;
sCritical = fmin(lambda x: -interpolate.splev(x, tckRetMap), 1)*0.998

def fCritical(s):
    po = retmapm(3, s) - s
    return po
sCritical = newton(fCritical, sCritical, tol=1.48e-12)


print "Scritical:"
print sCritical

Kneading = np.copy(sCritical)
KneadingSequence = ''
KneadingValueBin = '0.'
KneadingValue = 0
for i in range(nMax):
    xnext = retmap(Kneading[-1])
    Kneading=np.append(Kneading, xnext)
    if xnext > sCritical:
        KneadingSequence = KneadingSequence+'1'
        if i == 0:
            KneadingValueBin = KneadingValueBin+'1'
        else:
            KneadingValueBin = KneadingValueBin+str(int(not(int(KneadingValueBin[-1]))))
    else:
        KneadingSequence = KneadingSequence+'0'
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
    itinerary = ''
    for i in range(n):
        s = retmap(s)
        if s>sCritical:
            itinerary = itinerary + '1'
        elif s<sCritical:
            itinerary = itinerary + '0'
    return itinerary

def Splus2gamma(itinerary):
    gamma = 0
    for i in range(len(itinerary)):
        if i == 0:
            gammaBin = '0.'+itinerary[i]
        elif itinerary[i]=='0':
            gammaBin = gammaBin + gammaBin[-1]
        elif itinerary[i]=='1':
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

print "Kneading value fun:"
print TopologicalCoordinate(sCritical, 8)
    
if computeRPO:

    #Look upto the mth return map to find the periodic orbit candidates:
    scandidates = np.zeros([1,2]) #dummy matrix to hold po candidate arclengths
    smin = np.min(sn)
    smax = np.max(sn)
    
    for i in range(m):
        #Define the function zeros of which would correspond to the periodic 
        #orbit arclengths:
        def fpo(s):
            po = retmapm(i+1, s) - s
            return po
        #print "retmap %i" %(i+1)
        fpoevo=0
        for s0 in np.arange(smin, smax, (smax-smin)/50000):
            fpoev = fpo(s0)
            if fpoev * fpoevo < 0: #If there is a zero-crossing, look for the root: 
                sc = newton(fpo, s0, tol=1.48e-12)
                #print "sc = %f" %sc
                newcandidate = 1
                for j in range(np.size(scandidates,0)):
                    #Discard if found candidate is previously found:
                    if np.abs(scandidates[j,1] - sc)<1e-9: 
                        newcandidate=0
                if newcandidate:
                    scandidates = np.append(scandidates, np.array([[i+1, sc]], float), axis=0)
            fpoevo = fpoev
    
    scandidates = scandidates[1:np.size(scandidates,0), :]
    print "Periodic orbit candidate arclengths:"
    print scandidates
    
    #Read the maximum value of the return map to determine the point of the 
    #binary partition:
    imaxint = np.argmax(yintRetMap)
    sboundary = xintRetMap[imaxint]
    
    #Dummy itinerary array:
    scitineraries = np.array([''])
    #Determine itineraries:
    for i in range(np.size(scandidates,0)):
        scandidate = scandidates[i,1]
        if scandidate < sboundary:
            itinerary = "0"
        else:
            itinerary = "1"
        m = int(scandidates[i,0])   
        for i in range(1,m):
            scandidate = retmap(scandidate) 
            if scandidate < sboundary:
                itinerary = itinerary+"0"
            else:
                itinerary = itinerary+"1"
        scitineraries = np.append(scitineraries, np.array([itinerary]))
    scitineraries = scitineraries[1:]
    print "itineraries:"
    print scitineraries
    
    #Dummy assignments:
    group = np.zeros(len(scitineraries), int)
    position = np.zeros(len(scitineraries), int)
    #Group candidates into multiple shooting candidates according to their itineraries:
    for i in range(np.size(scandidates,0)):
        "ith itinerary:"
        iti = scitineraries[i]
        #print "range(length of iti - 1)"
        #print range(len(iti)-1)
        #Indices of mth order periodic orbit candidates:
        imth = np.argwhere(scandidates[:,0]==scandidates[i,0])
        if len(imth)==1:
            group[i] = int(1)
            position[i] = int(1)
        if group[i] == 0:
            group[i] = int(np.max(group))+1 #ith candidate starts a new group
            position[i] = int(1)
            for shift in range(-len(iti)+1,0):
                for j in imth:
                #print j
                    if scitineraries[j]  == (iti[shift:] + iti[0:len(iti)+shift]):
                        group[j] = group[i] #jth candidate belong to the same group
                        position[j] = shift + len(iti)+1                    
    print "group:"      
    print group     
    print "position:"
    print position
    together = np.append(np.array([scitineraries]).transpose(), np.array([group]).transpose(), axis=1)
    together = np.append(together, np.array([position]).transpose(), axis=1)
    together = np.append(together, np.array([scandidates[:,1]], str).transpose(), axis=1)
    print "Itinerary | Group | Position in the group"
    print together
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
        p0x = ps2D[i0,0]
        
        px = newton(fs2ps2D, p0x, args=((s*snmax)+snmin,), tol=1e-12)
        py = interpolate.splev(px, tckps)
        
        return np.array([px, py])
    
    ps2Dcandidates = np.array([s2ps2D(s) for s in scandidates[:,1]] )
    print "ps2Dcandidates:"
    print ps2Dcandidates
    
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
    
    pscandidates = np.array([ps2D2psxhat(ps2D) for ps2D in ps2Dcandidates] )
    print "pscandidates"
    print pscandidates
    
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
        print "Tapproximate:"
        print Tapproximate
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
        print "psreturn:"
        print psreturn
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
    
    TOF = np.array([timeofflight(x0) for x0 in pscandidates], float)
    print "time of flights:"
    print TOF
    phiret = np.array([phireturn(pscandidates[i,:], TOF[i]) for i in range(np.size(TOF,0))], float)
    print "phiret:"
    print phiret
    raw_input("Press Enter to continue...") 
    Adaptive = False
    iAdaptiveMax = 20           
    factor = 2
    
    tol = 1e-9
    earray = np.array([], float)
    xrpo = np.zeros(np.shape(pscandidates))
    tofrpo = np.zeros(np.shape(TOF))
    phirpo = np.zeros(np.shape(phiret))
    
    for i in range(1, int(np.max(group)+1)):
    
        AdaptiveFail = False
        print "Group no: ", i
        
        #Get corresponding indices for the group i:
        gindices = np.argwhere(group == i)
        #Reshape them in a 1D array:
        gindices = gindices.reshape(np.size(gindices))
        #Get coordinates of group i
        xi = pscandidates[gindices,:]
        #Get time of flights for group i:
        Ti = TOF[gindices]
        #Get group parameter for group i:
        phii = phiret[gindices]
        #Get the positions in group i
        posi = np.array([position[gindices]-1], int)
        posi = posi.reshape(np.size(posi))
        print "posi = ", posi 
        #Sort with respect to positions:
        xi = xi[posi, :]
        print "xi = ", xi 
        Ti = np.array([Ti[posi]])
        Ti = Ti.reshape(np.size(Ti))
        print "Ti = ", Ti 
        phii = np.array([phii[posi]])
        phii = phii.reshape(np.size(phii))%(2*np.pi)
        print "phii = ", phii 
        raw_input("Press Enter to continue...") 
        #How many points:
        npts = np.size(gindices)
        #Num of dim:
        n = 4
        #Initiate A, xT, xTT, error, error00 matrices:
        A = np.zeros(((n+2)*npts, (n+2)*npts))
        xT = np.zeros((npts, n))
    
        error = np.zeros((npts, n))
        error00 = np.zeros(npts*(n+2))
        
        #Dummy variables used in the adaptive block:
        xTT = np.zeros(np.shape(xT))
        xii = np.zeros(np.shape(xi))
        Tii = np.zeros(np.shape(Ti))
        phiii = np.zeros(np.shape(Ti))
        errorr = np.zeros(np.shape(error))
        
        #Define maximum number of iterations:
        itermax = 200
        #Apply ChaosBook p294 (13.11) 
        #with constraint nhat . Dx = 0
        iteration=0
        converged = False
        earray = np.array([], float)
        
        while not(converged):
            
            print "xi"
            print xi
            
            A = np.zeros(((n+2)*npts, (n+2)*npts))
            for j in range(npts):
                
                xTj = twomode.ftau(xi[j], Ti[j])
                xT[j, :] = xTj
                
                #Construct Aj matrix for each cycle element:
                #Jacobian:
                J = twomode.Jacobian(xi[j], Ti[j])
                
                #minusgvx term:
                vx = np.array(twomode.vfullssp(xTj,0), float)
                minusgvx = -np.dot(twomode.LieElement(phii[j]), vx)
                minusgvx = minusgvx.reshape(-1,1)
                
                #minusTgfTx term:
                minusTgfTx = np.dot(twomode.LieElement(phii[j]), xTj)
                minusTgfTx = -np.dot(T, minusTgfTx)
                minusTgfTx = minusTgfTx.reshape(-1,1)
                
                Aj = np.concatenate((-np.dot(twomode.LieElement(phii[j]), J), minusgvx), axis=1)
                Aj = np.concatenate((Aj, minusTgfTx), axis=1)
                
                #nhat0 = np.append(nhat.reshape(-1,1), np.array([[0], [0]], float))
                
                #print nhat0
                
                Aj = np.concatenate((Aj, np.append(vx.reshape(-1,1), np.array([[0], [0]], float), axis=0).transpose()), axis=0)
                #Aj = np.concatenate((Aj, np.append(nhat.reshape(-1,1), np.array([[0], [0]], float), axis=0).transpose()), axis=0)
                #Aj = np.concatenate((Aj, np.append(tp.reshape(-1,1), np.array([[0], [0]], float), axis=0).transpose()), axis=0)        
                tx = np.dot(T, xTj)
                Aj = np.concatenate((Aj, np.append(tx.reshape(-1,1), np.array([[0], [0]], float), axis=0).transpose()), axis=0)     
                
                #print Aj
    
                A[j*(n+2):(j+1)*(n+2), j*(n+2):(j+1)*(n+2)] = Aj
            
                A[j*(n+2):j*(n+2)+n, ((j+1)%npts)*(n+2) : ((j+1)%npts)*(n+2) + n] = A[j*(n+2):j*(n+2)+n, ((j+1)%npts)*(n+2) : ((j+1)%npts)*(n+2) + n] + np.identity(n)
            
            for k in range(npts):
                for l in range(npts):
                    print "A("+str(k)+", "+str(l)+")"
                    print A[k*(n+2):(k+1)*(n+2), l*(n+2):(l+1)*(n+2)]
                        
            for j in range(npts):
                
                error[j, :] = xi[j] - np.dot(twomode.LieElement(phii[(j-1)%npts]), xT[(j-1)%npts,:])
                error00[((j-1)%npts)*(n+2):((j-1)%npts + 1)*(n+2)] = np.append(error[j,:], np.array([0, 0]))
            
            print "error:"
            print error
            
            earray = np.append(earray, np.max(np.abs(error)))
                        
            Ainv = np.linalg.inv(A)
            
            print "Ainv.A"
            print np.dot(Ainv, A)
            #idntty = np.dot(Ainv, A)
        
            DxTphi = np.dot(Ainv, -error00)
            
            alpha = 1
            iadaptive = 0
            
            if Adaptive:
                
                Converging = False
                
                while not(Converging):
                    
                    DxTphi = DxTphi*alpha
                    iadaptive += 1
                    print "iadaptive= ", iadaptive
                    print "alpha= ", alpha
                    
                    for j in range(npts):
                
                        xii[j,:] = xi[j,:] + DxTphi[j*(n+2):j*(n+2)+n]
                        Tii[j] = Ti[j] + DxTphi[j*(n+2)+n]
                        phiii[j] = phii[j] + DxTphi[j*(n+2)+n+1]            
                    
                        xTT[j, :] = twomode.ftau(xii[j], Tii[j])
        
                    for j in range(npts):
                        
                        errorr[j, :] = xii[j] - np.dot(twomode.LieElement(phiii[(j-1)%npts]), xTT[(j-1)%npts,:])
                    
                    emax = np.max(np.abs(errorr))
                    emaxprevious = earray[len(earray)-1]
                        
                    if  emax < emaxprevious:
                        
                        Converging = True
                    
                    elif iadaptive > iAdaptiveMax:
                        
                        print "Adaptive step failed, looks kinda hopeless."
                        AdaptiveFail = True
                        raw_input("Press Enter to continue...") 
                        break
                    
                    else:                   
                        
                        alpha = float(alpha) / float(factor)
                        
                    
    
            for j in range(npts):
                #print "j ="    
                #print j    
                #print "((j-1)%npts)*(n+2) ="
                #print ((j-1)%npts)*(n+2)
                xi[j,:] = xi[j,:] + DxTphi[j*(n+2):j*(n+2)+n]
                #print "((j-1)%npts)*(n+2)+n ="
                #print ((j-1)%npts)*(n+2)+n
                Ti[j] = Ti[j] + DxTphi[j*(n+2)+n]
                #print "((j-1)%npts)*(n+2)+n+1 ="
                #print ((j-1)%npts)*(n+2)+n+1
                phii[j] = phii[j] + DxTphi[j*(n+2)+n+1]
    
    
            iteration += 1
            
            if np.max(np.abs(error)) < tol:
            
                xrpo[gindices, :] = xi
                tofrpo[gindices] = Ti
                phirpo[gindices] = phii
                converged = True
            
            if iteration == itermax or AdaptiveFail:
    
                if AdaptiveFail:
                    
                    print "Adaptive step failed, looks kinda hopeless."
                
                print "did not converged in given maximum number of steps"
                print "exitting..."
                xrpo[gindices, :] = xi
                tofrpo[gindices] = Ti
                phirpo[gindices] = phii
                break
            
        
        print "x_rpo:"
        print xrpo
        print "T_rpo:"
        print tofrpo
        print "phi_rpo:"
        print phirpo
        
        #plot(np.arange(iteration), earray)
        #xlabel('Iteration')
        #ylabel('Error')
        #savefig('image/errorrposearch'+str(i)+'.png', bbox_inches='tight', dpi=150)
        #show()
        
        #raw_input("Press Enter to continue...")

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
    xticks = np.linspace(min(ps2D[:,0]), max(ps2D[:,0]), Nticks)
    ax.set_xticks(xticks)
    ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=16); 
    yticks = np.linspace(min(ps2D[:,1]), max(ps2D[:,1]), Nticks)
    ax.set_yticks(yticks)
    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    yticks=yticks*100
    ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=16);  
    
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
    ax.set_aspect('equal')
    smin = np.min(sn)
    smax = np.max(sn)
    ax.set_xlim(smin,smax)
    ax.set_ylim(smin,smax)
    ax.set_xlabel('$s_n$', fontsize=24)
    ax.set_ylabel('$s_{n+1}$', fontsize=24)
    Nticks = 5

    xticks = np.linspace(smin, smax, Nticks)
    ax.set_xticks(xticks)
    ax.set_xticklabels(["$%.1f$" % xtik for xtik in xticks], fontsize=16); 

    yticks = np.linspace(smin, smax, Nticks)
    ax.set_yticks(yticks)
    ax.set_yticklabels(["$%.1f$" % ytik for ytik in yticks], fontsize=16); 
    
    #plt.figure(2, figsize=(8,8))
    #sp3 = np.array([retmapm(3, sn) for sn in srange])
    #plot(srange, sp3, 'b')
    #plt.hold(True)
    #plot(srange,srange,'g')
    
    savefig('RetMap.pdf', bbox_inches='tight', dpi=100) 
    call(["pdfcrop", "RetMap.pdf", "RetMap.pdf"], shell=True)
    
    show()
