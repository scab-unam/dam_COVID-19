# PDDE models
import numpy as np
from scipy.stats import binom
import matplotlib.pylab as gr

# -----------------------
# Numerics
# -----------------------
def RK2_Autonomous(f, p, parNames=[],parValues=[]):
    """Second-order Runge-Kutta method to solve x' = f(x) with U(t[0]) = U0.
    NOTES:
        This version is based on the algorithm presented in "Numerical
        Analysis", 6th Edition, by Burden and Faires, Brooks-Cole, 1997.
    """
    nForc=len(parNames)
    U=np.zeros((p['nSteps'], np.prod(np.shape(p['ic']))),"float64")
    U[0]=p['ic']
    if nForc>0:
        for i in range(p['nSteps']-1):
            for nn in range(nForc):
                p[parNames[nn]]=parValues[nn][i]
            k1 = p['stepSize'] * f( U[i], p) / 2.0
            U[i+1] = U[i] + p['stepSize'] * f( U[i] + k1, p)
    else:
        for i in range(p['nSteps']-1):
            k1 = p['stepSize'] * f( U[i], p) / 2.0
            U[i+1] = U[i] + p['stepSize'] * f( U[i] + k1, p)
    return U.transpose()

def sigmoid(u,u0=1.0/24.0,n=2):
    """ Usage:
    s = sigmoid(u,u0=1.0/24.0,n=2)
    """
    uu = u**n
    return uu/(uu+ u0**n)


def calcCFR(c,f):
    i = np.where(f==0)[0]
    ff = np.copy(f)
    ff[i]=1
    return c/ff

def geometricProb(p,n):
    return (1-p)**(n-1) * p


def dic2ClassVars(cla, dic):
    for k,v in dic.items():
        str="cla."+k+"= v"
        exec(str)
        print(str)
    return

def addSuffix(p,suff='',removeOriginal=1):
    pa= p.copy()
    for key, val in p.items():
        pa[key+suff]=val
        if removeOriginal:
            pa.pop(key)
    return pa

class rSECURD:
    def __init__(self, p):
        self.p = p
        print(self.p)
        self=dic2ClassVars(self,p)
        self.sampTimes= np.arange(0,p['timeMax'],p['timeStep'])
        self.nSteps = len(self.sampTimes)
        return

    def randSECURD(self,Z):
        S,E,C,U,R=Z;
        N = S+E+C+U+R
        pContact = (self.pExposure_E*E + self.pExposure_C* C +  self.pExposure_U*U)/N
        pTransmGivenExposure= self.pViremia1_contact * sigmoid(self.timeStep)
        susc = np.int32(S * self.pExposure_S)
        pE=pTransmGivenExposure * pContact
        newE = binom.rvs(susc, pE, size=1)
        newI= binom.rvs(E, self.pViremia2, size=1)
        newC = binom.rvs(newI, self.pConfirmed, size=1)
        newU = newI-newC
        deathC = binom.rvs(C,self.pDeathC,size=1)
        deathU = binom.rvs(U,self.pDeathU,size=1)
        newRC = binom.rvs(np.floor(self.pAccess*(C-deathC)),self.pHeal,size=1)
        newRU = binom.rvs(np.floor(self.pAccess*(U-deathU)),self.pHeal,size=1)
        SS = S - newE
        EE = E + newE - newI
        CC = C + newC - newRC - deathC
        UU = U + newU - newRU - deathU
        RR = R + newRC + newRU
        DDC = deathC
        DDU = deathU
        return [SS, EE, CC, UU, RR, DDC,DDU]

    def discreteDynamics(self, nonAutPars={}):
        """
        nAutPars is a dictionary containing parameters to be changed over time.
        The parameter names are the keys of the dictionary, the values are vectors with as many elements as the time series.
        """
        S=np.zeros(self.nSteps,'int32')
        E=np.zeros(self.nSteps,'int32')
        C=np.zeros(self.nSteps,'int32')
        U=np.zeros(self.nSteps,'int32')
        R=np.zeros(self.nSteps,'int32')
        DC=np.zeros(self.nSteps,'int32')
        DU=np.zeros(self.nSteps,'int32')
        print(self.nSteps)
        S[0] = self.S0; E[0]=self.E0
        C[0] = self.C0; U[0]=self.U0
        R[0] = self.R0; DC[0]=self.DC0; DU[0]=self.DU0
        for n in range(self.nSteps-1):
            for key,val in nonAutPars:
                self.p[key]=val[n]
            S[n+1],E[n+1],C[n+1],U[n+1],R[n+1],DC[n+1],DU[n+1] = self.randSECURD(Z=[S[n],E[n],C[n],U[n],R[n]])
        self.S = S; self.E=E;self.C=C;self.U=U;self.R=R;self.DC=DC; self.DU=DU
        return [S,E,C,U,R,DC,DU]

class dSECURD:
    def __init__(self, p):
        self.p = p
        self=dic2ClassVars(self,p)
        print(self)
        self.sampTimes= np.arange(0,self.p['timeMax'],self.p['timeStep'])
        self.nSteps = len(self.sampTimes)
        return

    def deteSECURD(self,Z):
        S,E,C,U,R,DC,DU=Z;
        N = S+E+C+U+R
        pContact = (self.pExposure_E*E + self.pExposure_C* C +  self.pExposure_U*U)/N
        pTransmGivenExposure= self.pViremia1_contact * sigmoid(self.timeStep)
        susc = np.int32(S * self.pExposure_S)
        pE=pTransmGivenExposure * pContact
        newE = susc*pE
        newI= E*self.pViremia2
        pRecov = self.pAccess**self.pHeal
        deathC = self.pDeathC * C
        deathU = self.pDeathU * U
        newRC = pRecov*(C-deathC)
        newRU = pRecov*(U-deathU)
        dS = - newE
        dE = newE - newI
        dC = newI*self.pConfirmed - newRC - deathC
        dU = newI*(1-self.pConfirmed) - newRU - deathU
        dR = newRC + newRU
        dDC = deathC
        dDU = deathU
        return [dS, dE, dC, dU, dR, dDC,dDU]

    def dynamics(self, nonAutPars={}):
        """
        nAutPars is a dictionary containing parameters to be changed over time.
        The parameter names are the keys of the dictionary, the values are vectors with as many elements as the time series.
        """
        S,E,C,U,R,DC,DU= RK2_Autonomous(f=self.deteSECURD, p=self.p, parNames=[],parValues=[])
        self.S = S; self.E=E;self.C=C;self.U=U;self.R=R;self.DC=DC; self.DU=DU
        return [S,E,C,U,R,DC,DU]

randomDynamics=1
if randomDynamics==1:
    ranP={'pExposure_S':0.5, 'pExposure_E':1.0, 'pExposure_C':0.1, 'pExposure_U': 0.5,
    'pViremia1_contact':0.3, 'pViremia2': 0.25,
    'pDeathU':0.0001,'pDeathC':0.05, 'pHeal':0.1, 'pAccess':1.0, 'pConfirmed': 1/50.0,
    'S0':3.3*(10**3)-1, 'E0':2,'C0':0,'U0':0,'R0':0,'DC0':0, 'DU0':0,
    'timeStep':1.0, 'timeMax':365.0}
    ranP['ic'] = np.array([ranP['S0'], ranP['E0'], ranP['C0'], ranP['U0'], ranP['R0'], ranP['DC0'], ranP['DU0']])
    #
    rEpi = rSECURD(p=ranP)
    #x = uh.dSECURD([2000,1000,0,0,0,0])
    S,E,C,U,R,DC,DU= rEpi.discreteDynamics()
    #
    rEpi.N = S+E+C+U+R
    rows=4; cols=1;
    f= gr.figure(figsize=(15,9)); gr.ioff()
    ax=list()
    for n in range( (rows * cols)):
        ax.append(f.add_subplot(rows,cols,n+1))
    #
    aSR,aECU,aD,aCFR = ax
    aSR.plot(np.float64(rEpi.S)/rEpi.N,'-',ms=2,alpha=1,label=r'$(t,S)$, proportion, N(0)=%d'%rEpi.N[0])
    aSR.plot(np.float64(rEpi.R)/rEpi.N,'-',ms=2,label=r'$(t,R)$ proportion')
    aECU.plot(rEpi.E,'-',ms=2,label=r'$(t,E)$')
    aECU.plot(rEpi.C,'-',ms=2,label=r'$(t,C)$')
    aECU.plot(rEpi.U,'-',ms=2,label=r'$(t,U)$')
    aD.plot(rEpi.DC,'-',ms=2,label=r'$(t,D_C)$ incidence')
    aD.plot(rEpi.DU,'-',ms=2,label=r'$(t,D_U)$ incidence')
    aCFR.plot(100*calcCFR((rEpi.DC+rEpi.DU).cumsum(),(rEpi.C+rEpi.U).cumsum()),'-',ms=2,label=r'% CFR total cases')
    aCFR.plot(100*calcCFR(rEpi.DC.cumsum(),rEpi.C.cumsum()),'-',ms=2,label=r'% CFR confirmed')
    for n in range(rows):
        ax[n].set_xlabel('time (days)')
        #ax[n].set_ylabel('# individuals');
        ax[n].set_xlim(0,rEpi.timeMax);
        ax[n].legend();
    gr.ion(); gr.draw(); gr.show()


detDynamics=1
if detDynamics==1:
    detP={'pExposure_S':0.5, 'pExposure_E':1.0, 'pExposure_C':0.1, 'pExposure_U': 0.5,
    'pViremia1_contact':0.3, 'pViremia2': 0.25,
    'pDeathU':0.0001,'pDeathC':0.02, 'pHeal':0.1, 'pAccess':1.0, 'pConfirmed': 1/50.0,
    'S0':3.3*(10**3)-1, 'E0':2,'C0':0,'U0':0,'R0':0,'DC0':0, 'DU0':0, 'timeStep':1/10.0, 'timeMax':365.0}
    detP['ic'] = np.array([detP['S0'], detP['E0'], detP['C0'], detP['U0'], detP['R0'], detP['DC0'], detP['DU0']])
    #
    dEpi = dSECURD(detP)
    #x = uh.dSECURD([2000,1000,0,0,0,0])
    S,E,C,U,R,DC,DU= dEpi.dynamics()
    #
    dEpi.N = S+E+C+U+R
    rows=4; cols=1;
    f= gr.figure(figsize=(15,9)); gr.ioff()
    ax=list()
    for n in range( (rows * cols)):
        ax.append(f.add_subplot(rows,cols,n+1))
    #
    aSR,aECU,aD,aCFR = ax
    aSR.plot(np.float64(dEpi.S)/dEpi.N,'-',ms=2,alpha=1,label=r'$(t,S)$, proportion, N(0)=%d'%dEpi.N[0])
    aSR.plot(np.float64(dEpi.R)/dEpi.N,'-',ms=2,label=r'$(t,R)$ proportion')
    aECU.plot(dEpi.E,'-',ms=2,label=r'$(t,E)$')
    aECU.plot(dEpi.C,'-',ms=2,label=r'$(t,C)$')
    aECU.plot(dEpi.U,'-',ms=2,label=r'$(t,U)$')
    aD.plot(dEpi.DC,'-',ms=2,label=r'$(t,D_C)$ incidence')
    aD.plot(dEpi.DU,'-',ms=2,label=r'$(t,D_U)$ incidence')
    aCFR.plot(100*calcCFR((dEpi.DC+dEpi.DU).cumsum(),(dEpi.C+dEpi.U).cumsum()),'-',ms=2,label=r'% CFR total cases')
    aCFR.plot(100*calcCFR(dEpi.DC.cumsum(),dEpi.C.cumsum()),'-',ms=2,label=r'% CFR confirmed')
    for n in range(rows):
        ax[n].set_xlabel('time (days)')
        #ax[n].set_ylabel('# individuals');
        ax[n].set_xlim(0,dEpi.timeMax);
        ax[n].legend();
    gr.ion(); gr.draw(); gr.show()
