# Mathematical epidemiology base Code
# Created MAHV, 202004170503

import numpy as np
from scipy.stats import binom,poisson
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

def expF(d,tau):
    return 1- np.exp(-d/tau)

# -----------------------
# -----------------------
# -----------------------
# General SIR-like dynamics
# -----------------------
# -----------------------
# -----------------------
def NIW(U,p):
    """
    N - Non infected
    I - Infectious
    W - Withdrawn
    """
    N,I,W=U
    T = N+I+W
    newInfections= p['beta']*N*I/T
    newWithdrawals= p['mu']*I
    dN = -newInfections
    dI = newInfections-newWithdrawals
    dW = newWithdrawals
    return dN, dI, dW

# -----------------------
# General NIW-like dynamics with 4 infectious groups and different levels of exposure depending on the group. The main idea is to quantify the overall number of cases regardless of their condition. For questions related to use of resources, dynamics between the different stages of infection should be proposed.
# -----------------------
def pContact(exposure,sizes):
    return np.dot(exposure,sizes)

def randNIWIncidence(Z,p):
    """
    Motivated by the COVID-19 dynamics. Takes into account the infectious period of different risk groups according to their typical clinical histories.
    """
    #print(len(Z),len(Z[2]))
    N,W,I = Z;
    #print(I)
    INon,IMild,ISevere,IFatal = I
    T = N+I.sum()+W
    pExposure = np.dot( p['exposure_I'], I )/ T
    #print(pExposure)
    susc = N * p['exposure_N'] * p['weight_I']
    newI = np.zeros(p['nI'],'int64')
    newW_I=np.zeros(p['nI'],'int64')
    if susc.min()<10**3:
        for nn in range(p['nI']):
            newI[nn]= binom.rvs( np.int32( susc[nn] ), p['tGe'] * pExposure,loc=0, size=1)
            newW_I[nn]= binom.rvs(I[nn], p['pWithdraw'][nn],loc=0, size=1)
    else:
        for nn in range(p['nI']):
            newI[nn]= poisson.rvs( susc[nn] * p['tGe'] * pExposure, size=1)
            newW_I[nn]= poisson.rvs(I[nn] * p['pWithdraw'][nn], size=1)

    NN = - newI.sum()
    II = + newI - newW_I # n-dimensional vector
    WW = newW_I.sum()
    return np.array([NN,WW, II])


def randNIWDynamics(p, rhs=randNIWIncidence, nonAutPars={}):
    """
    nAutPars is a dictionary containing parameters to be changed over time.
    The parameter names are the keys of the dictionary, the values are vectors with as many elements as the time series.
    """
    Z = list()
    Z.append(np.array(p['iC']))
    for n in range(p['nSteps']-1):
        for key,val in nonAutPars:
            p[key]=val[n]
        Z.append(Z[n] + rhs(Z[n],p))
    return np.array(Z).transpose()





def uvField(ax,beta=0.5, delta=0.1, parts=100j):
    V,U=sc.mgrid[0:1:parts, 0:1:parts]
    dw=delta*V
    du= -beta*U*V
    dv= -du - dw
    #speed = sc.sqrt( U*U + V*V)
    speed = 2*sc.sqrt(V*V)
    #lw=1;
    db=delta/beta
    print('delta / beta = %f'%(db))
    lw=2*speed/speed.max()
    ax.plot([db,db],[0,1],'k--',lw=1)
    #ax.streamplot(U,V, du,dv, color='k', linewidth=lw, cmap=gr.cm.gray_r)
    #ax.streamplot(U,V, u,v, color=v, linewidth=2, cmap=gr.cm.autumn)
    #ax.streamplot(U,V, u,v, color='k', density=[0.75, 0.75], linewidth=lw)
    ax.streamplot(U,V, du,dv, color='k', linewidth=lw)
    return V,U,du,dv
