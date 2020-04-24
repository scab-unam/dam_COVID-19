#!/usr/bin/env python
# coding: utf-8

# # Probabilistic epidemiological model based on binomial sampling
# ### Carlos Ignacio Herrera-Nolasco$^{1}$ & Marco Arieli Herrera-Valdez$^{1}$
# #### $^{1}$ Laboratorio de Fisiología de Sistemas, Departamento de Matemáticas, Facultad de Ciencias, UNAM
#
# Consider a population of size $N(t)$ at time $t$.
# To model infection dynamics with a probabilistic approach, let $S$ and $I$, represent the number of susceptible and  infected individuals, respectively. Also, let $R$ represent those individuals who do not participate in the chain of transmission, either because they are immune, or because they died of infection. Then $N=S+I+R$ at each point in time. Let $x(t)$ is a random variable that represents the number of newly infected people between times $t$ and $t+\delta$, where $\delta$ is a small time step. The random variable $x$ should depend, at each time $t$, on the number of people susceptible to infection, and their probability of becoming infected. Similarly, let $y(t)$ be a random variable representing the number of people removed from the chain of infection at time $t$, which should depend on the current number of infected people $I(t)$ and the probability $\gamma$ of being removed from the chain of infection.
# The dynamics for $S$ and $I$ can then be written as
# \begin{eqnarray}
# S(t+\delta) &=& S(t) -  x(t),
# \\
# I(t+\delta)  &=& I(t) + x(t) - y(t).
# \end{eqnarray}
#
# For instance, if susceptible people are assumed to become infected when having an infectious encounter independently of one another, then $x(t)$ can be assumed to be a binomial random variable with the number of trials equal to $S(t)$. The probability of infection per unit time can then be expressed as the product of the probability of having contact with an infectious individual times the probability of successful transmission of the infection. If homogeneous mix is assumed, then $P(\textrm{infectious contact per unit time})=I(t)/N$. However, in a more realistic setting where mixing in the population is not homogeneous, it should be the case that  $ P(\textrm{infectious contact per unit time})=\epsilon(t) I(t)/N$ where $\varepsilon(t) \in [0,1]$ is an adjustment factor representing the proportion of infected that susceptibles are exposed to at time $t$. As a result, $x(t) \sim Bin(S(t),\delta~\alpha(t))$ with
# \begin{equation}
# \alpha(t) = P(\textrm{infection at time } t) = \varepsilon(t) P(\textrm{infectious contact per unit time}) P(\textrm{successful transmission per contact}).
# \end{equation}
# where $\varepsilon$ the probability of exposure.
#
# Similarly, the infected individuals can be assumed to become removed  from the chain of transmission independently, which means that $y(t)\sim Bin(S(t),\delta~\gamma(t))$.

import scipy as sc
import numpy as np
from scipy.stats import binom
import matplotlib.pylab as gr

# ----------------------------------------------------------------
# Random dynamics for SIR with simple binomial sampling and exponential trigger times
# ----------------------------------------------------------------
class pSIR:
    def __init__(self, p):
        self.p = p
        self.dic2ClassVars(p)
        self.sampTimes= np.arange(0,p['timeMax'],p['timeStep'])
        self.nSteps = len(self.sampTimes)
        return

    def prepNumerics(self):
        self.sampTimes= np.arange(0,p['timeMax'],p['timeStep'])
        self.nSteps = len(self.sampTimes)
        self.ic = np.array([self.S0,self.I0,self.R0])
        return

    def dic2ClassVars(self, dic):
        for k,v in dic.items():
            str="self."+k+"= v"
            exec(str)
        return

    def pSIR(self,p,U):
        S,I,R=U
        N= S+I+R
        eS = sc.int32(S* p['exposure_S'])
        eI = sc.int32(p['exposure_I']*I)
        aI = sc.int32(p['accessHC_I']*I)
        contactEffectiveness=(1-np.exp(-p['timeStep']/p['infecContactTime']))
        pContact = contactEffectiveness * eI/N
        pInfec = p['timeStep']*p['pSuccTrans']*pContact
        pDeath = (1-np.exp(-p['timeStep']/p['deadTime']))
        pRemov = 1-np.exp(-p['timeStep']/p['removTime'])
        newI = binom.rvs(eS, pInfec, size=1)
        newD = binom.rvs(I, pDeath, size=1)
        newR = binom.rvs(eI, pRemov, size=1)
        dS = - newI
        dI = newI - newR
        dR = newR
        return np.array([dS,dI,dR])

    def discreteDynamics(self,p,parNames=[],parVals=[]):
        nForc=len(parNames)
        Z=np.zeros((p['nSteps'], np.prod(np.shape(self.ic))),"float64")
        Z[0]= self.ic
        for n in range(p['nSteps']-1):
            if nForc>0:
                for nn in range(nForc):
                    p[parNames[nn]]=parValues[nn][i]
            oo = self.pSIR(p,Z[n])
            Z[n+1] = Z[n] + oo[:,0]
        return Z.transpose()

p= {'nSteps':365, 'timeStep':1.0, 'timeMax':10, 'N':1000, 'S0':1000, 'I0':1, 'R0':0,
'infecContactTime':1.0/3600, 'exposure_S':1.0, 'exposure_I':1.0,'accessHC_I':1.0, 'pSuccTrans':0.5,
'deadTime':20.0,'removTime':7.0,'pDead':0.02}

z = pSIR(p)
z.prepNumerics()
if 0:
    S,I,R= z.pSIR(z.p,z.ic)
    S,I,R=z.discreteDynamics(p)


rows=2; cols=2
f= gr.figure(figsize=(15,7)); gr.ioff()
axS = f.add_subplot(rows,cols,1)
axI = f.add_subplot(rows,cols,3)
axSI = f.add_subplot(1,cols,2)
axS.plot(S,'.',ms=2,alpha=1,label=r'$(t,S)$')
axI.plot(I,'.',ms=2,label=r'$(t,I)$')
axSI.plot(S,I,'.',lw=1,alpha=1
          ,label=r'$(S,I)$')
#gr.plot(R,'.',ms=2,label=r'$R$')
axS.plot([0,p['nSteps']],[0,0],'k',lw=3,alpha=0.2)
axI.plot([0,p['nSteps']],[0,0],'k',lw=3,alpha=0.2)
axI.set_xlabel('time (days)')
axS.set_ylabel('Population size'); axI.set_ylabel('Population size')
axS.set_xlim(0,p['tMax']);axI.set_xlim(0,p['tMax'])
axS.legend(); axI.legend(); axSI.legend()
gr.ion(); gr.draw(); gr.show()
print(np.array([S,I,R]).transpose())

# ## Behaviour of the model for different population sizes
#
# A qualitative comparison of the dynamics for different population sizes shows that the dynamics tend to behave more deterministically as $N$ grows.

# In[130]:


lS = list(); lI=list(); lR=list()

p= {'nSteps':200,'timeStep':1.0,'exposure':0.5,'pSuccTrans':0.5,'pRemov':0.1, 'N':100}
p['S0/N']= 0.999
p['R0/N']= 0.0
popSizes= [100,1000,10000,100000]
for N in popSizes:
    p['N']=N
    p['S0']= sc.floor(p['N']*p['S0/N'])
    p['R0']= sc.floor(p['N']*p['R0/N'])
    p['I0']= p['N']-p['S0']-p['R0']
    print('Initial conditions: (S0,I0,R0)=(%d,%d,%d)'%(p['S0'],p['I0'],p['R0']))
    [S,I,R]=pSIR(p)
    lS.append(S); lI.append(I)
    print('State after %d steps: (S,I,R)=(%d,%d,%d)'%(p['nSteps'],S[-1],I[-1],R[-1]))

rows=2; cols=2
f= gr.figure(figsize=(17,11)); gr.ioff()
axS = f.add_subplot(rows,cols,1)
axI = f.add_subplot(rows,cols,3)
axSI = f.add_subplot(1,cols,2)
for n in range(len(popSizes)):
    axS.plot(lS[n]/sc.float32(popSizes[n]),'-',ms=2,alpha=0.65,label=r'$S_{%d}$'%popSizes[n])
    axI.plot(lI[n]/sc.float32(popSizes[n]),'-',ms=2*n,alpha=0.65,label=r'$I_{%d}$'%popSizes[n])
    axSI.plot(lS[n]/sc.float32(popSizes[n]),lI[n]/sc.float32(popSizes[n]),'-',ms=2*n,alpha=1,label=r'$(S_{%d},I_{%d})$'%(popSizes[n],popSizes[n]))
    axI.set_xlabel('time (days)')
    axS.set_ylabel('Normalized population size')
    axI.set_ylabel('Normalized population size')
    axSI.set_xlabel('Normalized population size')
    axS.legend(); axI.legend(); axSI.legend()
axS.plot([0,p['nSteps']],[0,0],'k',lw=1,alpha=0.2)
axI.plot([0,p['nSteps']],[0,0],'k',lw=1,alpha=0.2)
gr.ion(); gr.draw()
#print(sc.array([S,I,R]).transpose())


# ## Dynamics in small populations

# #### Variations due to randomness assuming that half of the susceptible population is exposed to the infection per unit time

# In[68]:


lS = list(); lI=list(); lR=list()
p= {'nSteps':200,'timeStep':1.0,'exposure':0.5,'pSuccTrans':0.5,'pRemov':0.1,'N':100}
p['S0/N']= 0.95
p['R0/N']= 0.0
p['S0']= sc.floor(p['N']*p['S0/N'])
p['R0']= sc.floor(p['N']*p['R0/N'])
p['I0']= p['N']-p['S0']-p['R0']
nRuns=5
print('Initial conditions: (S0,I0,R0)=(%d,%d,%d)'%(p['S0'],p['I0'],p['R0']))
for xx in range(nRuns):
    [S,I,R]=pSIR(p)
    lS.append(S); lI.append(I)
    print('State after %d steps: (S,I,R)=(%d,%d,%d)'%(p['nSteps'],S[-1],I[-1],R[-1]))
    print('Peak # of infected',I.max())
    print('-----------------------')

rows=2; cols=1
f= gr.figure(figsize=(15,7)); gr.ioff()
axS = f.add_subplot(rows,cols,1)
axI = f.add_subplot(rows,cols,2)
for n in range(nRuns):
    axS.plot(lS[n],'-',ms=2,alpha=0.65,label=r'$S_{%d}$'%n)
    axI.plot(lI[n],'-',ms=2*n,alpha=0.65,label=r'$I_{%d}$'%n)
    axI.set_xlabel('time (days)')
    axS.set_ylabel('Normalized population size')
    axI.set_ylabel('Normalized population size')
    axS.legend(); axI.legend()
axS.plot([0,p['nSteps']],[0,0],'k',lw=1,alpha=0.2)
axI.plot([0,p['nSteps']],[0,0],'k',lw=1,alpha=0.2)
gr.ion(); gr.draw()
#print(sc.array([S,I,R]).transpose())


# #### Variation due to randomness assuming that 7/10 of the susceptible population is exposed to the infection per unit time

# In[70]:


lS = list(); lI=list(); lR=list()
p= {'nSteps':200,'timeStep':1.0,'exposure':0.7,'pSuccTrans':0.5,'pRemov':0.1,'N':100}
p['S0/N']= 0.95
p['R0/N']= 0.0
p['S0']= sc.floor(p['N']*p['S0/N'])
p['R0']= sc.floor(p['N']*p['R0/N'])
p['I0']= p['N']-p['S0']-p['R0']
nRuns=5
print('Initial conditions: (S0,I0,R0)=(%d,%d,%d)'%(p['S0'],p['I0'],p['R0']))
for xx in range(nRuns):
    [S,I,R]=pSIR(p)
    lS.append(S); lI.append(I)
    print('State after %d steps: (S,I,R)=(%d,%d,%d)'%(p['nSteps'],S[-1],I[-1],R[-1]))
    print('Peak # of infected',I.max())
    print('-----------------------')

rows=2; cols=1
f= gr.figure(figsize=(15,7)); gr.ioff()
axS = f.add_subplot(rows,cols,1)
axI = f.add_subplot(rows,cols,2)
for n in range(nRuns):
    axS.plot(lS[n],'-',ms=2,alpha=0.65,label=r'$S_{%d}$'%n)
    axI.plot(lI[n],'-',ms=2*n,alpha=0.65,label=r'$I_{%d}$'%n)
    axI.set_xlabel('time (days)')
    axS.set_ylabel('Normalized population size')
    axI.set_ylabel('Normalized population size')
    axS.legend(); axI.legend()
axS.plot([0,p['nSteps']],[0,0],'k',lw=1,alpha=0.2)
axI.plot([0,p['nSteps']],[0,0],'k',lw=1,alpha=0.2)
gr.ion(); gr.draw()
#print(sc.array([S,I,R]).transpose())


# #### Variation due to randomness assuming that half of the susceptible population is exposed to the infection per unit time

# In[136]:


# Probabilistic SIR with simple binomial dynamics for the newly infected and the newly removed from the chain of infection
def pSIR_new(p):
    S=sc.zeros(p['nSteps'],'int32')
    I=sc.zeros(p['nSteps'],'int32')
    R=sc.zeros(p['nSteps'],'int32')
    newI=sc.zeros(p['nSteps'],'int32')
    newR=sc.zeros(p['nSteps'],'int32')
    S[0]=p['S0']; I[0]=p['I0'];R[0]=p['R0']
    N0=p['S0']+p['I0']+p['R0']

    for n in range(p['nSteps']-1):
        pContact = p['timeStep']*p['exposure']*I[n]/sc.float32(N0)
        pInfec = p['timeStep']*p['pSuccTrans']*pContact
        newI[n] = binom.rvs(S[n],pInfec,size=1)
        newR[n] = binom.rvs(I[n],p['pRemov'],size=1)
        #newR[n] = sc.int32(p['gamma']*I[n])
        S[n+1]=S[n] - newI[n]
        I[n+1]=I[n] + newI[n]-newR[n]
        R[n+1]=R[n] + newR[n]
    return [S,I,R,newI,newR]


# ### Total number of infected, peak time for infection, and maximum infected during the epidemic for small populations.
# It is assumed that half of the susceptible population is exposed to infection per unit time

# In[89]:


lS = list(); lI=list(); lR=list(); lnewI=list(); lnewR=list()
p= {'nSteps':200,'timeStep':1.0,'exposure':0.5,'pSuccTrans':0.5,'pRemov':0.1,'N':100}
p['S0/N']= 0.99
p['R0/N']= 0.0
p['S0']= sc.floor(p['N']*p['S0/N'])
p['R0']= sc.floor(p['N']*p['R0/N'])
p['I0']= p['N']-p['S0']-p['R0']
nRuns=5
print('Initial conditions: (S0,I0,R0)=(%d,%d,%d)'%(p['S0'],p['I0'],p['R0']))
for xx in range(nRuns):
    [S,I,R,newI,newR]=pSIR_new(p)
    lS.append(S); lI.append(I); lnewI.append(newI); lnewR.append(newR)
    print('State after %d steps: (S,I,R)=(%d,%d,%d)'%(p['nSteps'],S[-1],I[-1],R[-1]))
    print('Imax = %d, totalI=%d, peak time=%d days'%(I.max(),newI.sum(),I.argmax()))
    print('-----------------------')

rows=1; cols=1
f= gr.figure(figsize=(15,7)); gr.ioff()
#axS = f.add_subplot(rows,cols,1)
axI = f.add_subplot(rows,cols,1)
for n in range(nRuns):
    axI.plot(lI[n],'-',ms=2*n,alpha=0.8,label=r'$totalI=%d, I_{peak}=%d, t_{peak}=%d$ days'%(lnewI[n].sum(),lI[n].max(),lI[n].argmax()))
    axI.set_xlabel('time (days)')
    axS.set_ylabel('Normalized population size')
    axI.set_ylabel('Normalized population size')
    axS.legend(); axI.legend()
axS.plot([0,p['nSteps']],[0,0],'k',lw=1,alpha=0.2)
axI.plot([0,p['nSteps']],[0,0],'k',lw=1,alpha=0.2)
gr.ion(); gr.draw()
#print(sc.array([S,I,R]).transpose())


# In[ ]:
