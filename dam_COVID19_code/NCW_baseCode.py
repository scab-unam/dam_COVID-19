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




def plotEpidemic(N,W,I,p,casesData, maxCases=10000,maxDeaths=1000):
    INon= np.array([I[m][0] for m in range(len(I))])
    IMild= np.array([I[m][1] for m in range(len(I))])
    ISevere= np.array([I[m][2] for m in range(len(I))])
    IFatal= np.array([I[m][3] for m in range(len(I))])
    ITot = INon+IMild+ISevere+IFatal
    dead = list()
    for nn in range(p['nI']):
        dead.append(p['pDeath'][nn]*W)
    days = np.arange(0,len(N)) +p['offset']
    daysMexico = np.arange(0,len(MexicoCases_April16))
    peakInd =ITot.argmax()
    peakDay = days[peakInd]

    fig= gr.figure(figsize=(11,5)); gr.ioff()
    rows=1;cols=1; ax=list()

    for rc in range(rows*cols):
        ax.append(fig.add_subplot(rows,cols,rc+1))

    cax=inset_axes(parent_axes=ax[0],
                            width="30%", # width = 30% of parent_bbox
                            height="40%", # height : 1 inch
                            loc='center right')
    #dax=inset_axes(parent_axes=ax[2],
    #                        width="30%", # width = 30% of parent_bbox
    #                        height="50%", # height : 1 inch
    #                        loc='center left')
    #ax[1].plot(days, N,label='$N$ (no infectados)')
    #ax[1].plot(days, W,label='$W$ (recuperados)')
    ax[0].plot(days, INon,label=r'$I_{Non}$')
    ax[0].plot(days, IMild,label=r'$I_{Mild}$')
    ax[0].plot(days, ISevere,label=r'$I_{Severe}$')
    ax[0].plot(days, IFatal,label=r'$I_{Fatal}$')
    ax[0].plot(days, ITot,label=r'$Tot Cases$')
    str1=r'Confirmed cases x %d, S. Salud Mexico (04/16/2020)'%p['underReport']
    ax[0].plot(daysMexico, p['underReport']*casesData,'o',ms=2,label=str1)
    cax.plot(days, INon)
    cax.plot(days, IMild)
    cax.plot(days, ISevere)
    cax.plot(days, IFatal)
    cax.plot(days, ITot)
    cax.plot(daysMexico, p['underReport']* MexicoCases_April16,'o',ms=2,)
    #for rc in range(p['nI']):
    #    ax[2].plot(days, dead[rc],label=r'$Fallecidos_{%s}$'%p['stagesISpanish'][rc])
    #ax[2].plot(days, dead[rc].sum(0),label=r'$Fallecidos_{%s}$'%p['stagesISpanish'][rc])
    #dax.plot(days, dead[rc].sum(0))
    #ax[2].plot(daysMexico, underReport*MexicoDeaths_April16,label=r'$Fallecidos_{04/16}$')
    #dax.plot(daysMexico, underReport*MexicoDeaths_April16,label=r'$Fallecidos_{04/16}$')
    for rc in range(rows*cols):
        ax[rc].legend()
    cax.set_xlim(0,50); cax.set_ylim(0,maxCases)
    ax[0].set_ylabel('Confirmed cases')
    ax[0].set_xlabel('Days from first report on Feb 27, 2020')
    #dax.set_xlim(0,50); dax.set_ylim(0,maxDeaths)
    strTit= '''Covid-19 dynamics from the contribution of infectious individuals with different contagion intervals \n
    (Herrera-Nolasco, Herrera-Valdez, 2020)'''
    ax[0].text(peakInd+5, ITot.max()*.95, '%d days, %d infected'%(peakDay,ITot[peakInd]))
    ax[0].text(peakInd+5, ITot.max()*.9, 'Peak estimated around May 13, 2020')
    fig.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.85,wspace=0.2,hspace=0.25)
    fig.suptitle(strTit)
    gr.ion(); gr.draw(); gr.show()
    return fig


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
