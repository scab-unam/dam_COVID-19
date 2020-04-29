#!/usr/bin/env python
# coding: utf-8

# ## Analysis of pandemic dynamics assuming SIR-like evolution
# ### Carlos Ignacio Herrera-Nolasco$^{1}$, Marco Arieli Herrera-Valdez$^{1}$
# #### $^{1}$ Laboratorio de Fisiología de Sistemas, Facultad de Ciencias, UNAM
# 
# 
# ### Epidemic dynamics 
# 
# Assume that the population is divided in three subsets representing the non-infected and susceptible, the infected, and those that can no-longer participate in the chain of infections due to immunity or death. 
# Let the densities of those subpopulations be represented by $x$, $y$, and $z$ respectively, with $1=x+y+z$. As a consequence,  $\partial_t x = -\partial_t y - \partial_t z$. Let  $\alpha$ represent the infection rate given an infectious contact, and $\tau$ the average waiting time until an individual is no longer infected. The dynamics for $y$ can be written as 
# \begin{equation}
# \partial_t y = \alpha y \bigg( 1- z - y \bigg) - \beta y,
# \end{equation}
# and the evolution for $z$ is then
# \begin{equation}
# \partial_t z = \beta y
# \end{equation}
# 
# The incidence of cases is the number of new infected. From the equation for the change in $y$, the incidence is given by  $$f(y,z;\alpha) =\alpha y \bigg( 1- z - y \bigg),$$ which means that the incidence in the SIR model is a cuadratic function of the prevalence, but it also depends on $z$, the density of no-longer infected people. Since $z$ is an increasing function, then the graph of the incidence can be thought of as a curve in 3D space with quadratic shape. 
# 
# 
# #### Parameters for modeling 
# If $\gamma$ is the rate of removal by recovery and acquisition of immunity and $\delta$ is the fatality rate due to infection, we can try to obtain the parameters for the model from the data. To do so, consider the equation for $z$ to obtain,
# \begin{equation}\beta = \frac{\partial_t z}{y},\end{equation} 
# and also 
# \begin{equation}\alpha = \frac{\partial_t y + \partial_t z}{y \bigg( 1 -z - y  \bigg)}. \end{equation} 
# The problem is then to link the data to the variables $y$ and $z$.

# ### Linking the data to the modeling variables 
# 
# If $h$ is the time step for sampling, let $Y_{n}$ represent the number of cases at time $t=nh$. Note that $Y_n$ may be thought of in terms of the density of cases $y_n$. Explicitly, $y_n = Y_n/T $ here $T$ is the size of the local population.  From data, the cumulative number of cases at step $nh$ can be written as 
# $$C_{n} = \sum_{k=0}^{n} Y_k$$
# which means that the number of new cases (incidence) at time nh is 
# $$Y_n = C_n - C_{n-1}.$$
# The cumulative recoveries and deaths are 
# $$Z_n = R_{n}+D_{n}.$$
# The incidence is then 
# $$C_n - R_n - D_n$$
# 
# #### Parameters for modeling 
# We can try to obtain the parameters for the model from the data by approximated by substitution of the discrete differences 
# $$
# \partial_t y \approx \frac{u(t+h)-u(t)}{h},  \quad u \in \left\{Y,Z \right\}
# $$
# The parameters $\alpha$ and $\beta$ can then be estimated from the data. 
# 

# In[6]:


#import sys
#sys.path.insert(1, './')
#print(sys.path)
from dam_COVID19_baseCode import *
import matplotlib.pylab as gr
small={'family' : 'normal','weight' : 'normal','size'   : 8}
medium={'family' : 'normal','weight' : 'normal','size'   : 10}
large={'family' : 'normal','weight' : 'bold','size'   : 13}
gr.rc('font', size=small['size'], weight='normal')          # controls default text sizes
gr.rc('axes', titlesize=medium['size'])     # fontsize of the axes title
gr.rc('axes', labelsize=medium['size'])    # fontsize of the x and y labels
gr.rc('xtick', labelsize=small['size'])    # fontsize of the tick labels
gr.rc('ytick', labelsize=small['size'])    # fontsize of the tick labels
gr.rc('legend', fontsize=small['size'])    # legend fontsize
gr.rc('figure', titlesize=large['size'])  # fontsize of the figure title
get_ipython().magic(u'matplotlib inline')


# In[7]:


cases, deathCases,recovCases = getCSSEGISandData(urlData=1)


# In[8]:


cases.head(10)


# In[493]:


# ------------------------------
# Description of the data so that the headers and the columns
# without case data are distinguished
# ------------------------------
nRows,nCols=cases.shape
cases.head(10)
nHeaderRows=1;
#
nHeaderCols=3
# how to generate date lists from a baseline using the datetime
dates = cases.columns[4:]
nDays = len(dates)
days = np.arange(nDays)
print('Got data from %d days between %s and %s'%(nDays,dates[0],dates[-1]))
# -------------------
print("""""")
# -------------------
npCases = cases.to_numpy()
countries = np.unique(npCases[:,1])
nCountries = len(countries)
print('Considering data from {d} countries'.format(d=nCountries))
# -------------------
# Sum the counts from each country and construct a new array
# -------------------
# These arrays have the same size as the countries array (unique countries)
totCases=gatherDataByCountry(df=cases,nHeaderCols=4)
totDeathCases=gatherDataByCountry(df=deathCases,nHeaderCols=4)
totRecovCases=gatherDataByCountry(df=recovCases,nHeaderCols=4)
# Save all into a dictionary
G = {'cCases':totCases.transpose(), 'cDeaths':totDeathCases.transpose(),'cRecovs':totRecovCases.transpose(), 'countries':countries.transpose()}
# -------------------
# Search regions to illustrate the case-fatality ratios
# -------------------
Pops_Millions = {'Algeria':43851044, 'Argentina':45195774,'Australia':25499884,'Belgium':11433256,
                 'Brazil':212559417, 
                 'Canada':37742154, 'China':1439323776, 'Colombia':50882891,'Egypt':102334404,
                 'France':67886011,'Germany':83783942, 'Japan':126476461,'Korea, South':51269185,
                 'Indonesia':273523615, 'Iran':83992949,'Israel':8655535,'Italy':60461826,
                 'Mexico':128932753, 'Niger':24206644, 'Singapore':5850342, 
                 'South Africa':59308690,'Spain':46754778,
                 'United Kingdom':67886011, 'US':331002651, 'Venezuela':28870195}


# In[494]:


print(type(totCases))
print(totCases.shape)


# #### Parabola describing the incidence as a function of the cases

# In[495]:


def getSingleCountryData(G, popSize, cou,b = 0.002):
    cInd = np.where(G['countries']==cou)[0][0]
    cases= np.float64(G['cCases'][cInd,:])
    recovs= np.float64(G['cRecovs'][cInd,:])
    deaths= np.float64(G['cDeaths'][cInd,:])
    print(deaths.shape)
    country= {'cCases':cases,'cDeaths':deaths,'cRecovs':recovs, 'Z':deaths+recovs}
    country.update({'popSize':popSize[cou], 'cdCases':cases/popSize[cou], 'cdDeaths':deaths/popSize[cou], 'cdRecovs':recovs/popSize[cou], 'name':cou})
    return country
    
def dataSIR_Country(country,a=80,b=0.0002):
    # Prevalence time series
    # The time series of cumulative cases without cumulative recovs and without deaths
    Y = country['cCases'] - country['cRecovs'] - country['cDeaths'] 
    y = country['cdCases'] - country['cdRecovs'] - country['cdDeaths'] 
    # "Prevalence" of deaths and recoveries 
    z = country['cdRecovs'] + country['cdDeaths']
    Z = country['cRecovs'] + country['cDeaths']
    dY=np.zeros(len(Y)); dY[1:]= Y[1:]-Y[:1]
    dZ=np.zeros(len(Z)); dZ[1:]= Z[1:]-Z[:1]
    dy=np.zeros(len(y)); dy[1:]= y[1:]-y[:1]
    dz=np.zeros(len(z)); dz[1:]= z[1:]-z[:1]
    incidence = dY-dZ
    x = 1-z-y
    country.update({'x':x,'y':y,'z':z,'preval':Y, 'Z':Z,'dy':dy,'dz':dz,'dY':dY,'dZ':dZ,'incid':incidence})
    country['beta']= dz/y
    country['alpha']=(dy+dz)/(x*y)
    country['R_t']= country['alpha']*x/country['beta']
    country['fy'] = a*y*(b -y) 
    country['fY'] = country['popSize']*country['fy'] 
    country['Deltay']= a*y*(1-z -b -y) 
    country['DeltaY']= country['popSize']* a*y*(1-b-z -y) 
    return country

def beta_Estimate(country):
    country['beta']= country['dZ']/country['preval']
    return country

def alpha_Estimate(country):
    country['alpha']=(country['dY']+country['dZ'])/((1-country['Z']-country['preval'])*country['preval'])
    return country

def incidenceFit(country,a=1,b=0.0002):
    # Incidence 
    country['fy'] = a*country['y']*(b-country['y']) 
    country['fY'] = country['popSize']*country['fy'] 
    return country


# In[498]:


def plotDataPhasePortraits(country):
    maxY=1.1*country['preval'].max()
    ff =gr.figure(figsize=(17,5)); gr.ioff(); rows=1; cols=3
    ax1= ff.add_subplot(rows,cols,1); tax1= ax1.twinx()
    ax2= ff.add_subplot(rows,cols,2); #tax2= ax2.twiny()
    ax3= ff.add_subplot(rows,cols,3); tax3= ax3.twinx()
    tax1.plot(days, country['R_t'],'o',markeredgecolor='blue', markerfacecolor='white',label=r'$y_k$')
    tax1.plot([days[0],days[-1]], [1,1],'k:')
    ax1.plot(days, country['preval'],'b.',label=r'$R_t =\alpha (1-z-y)/\beta$')
    ax1.plot(days, country['cDeaths'],'.',label=r'$D_k$')
    ax1.plot(days, country['cRecovs'],'.',label=r'$R_k$')
    ax1.plot(days, country['Z'],'.',label=r'$Z_k=D_k + R_k$')
    ax2.plot(country['fy'],country['preval'],'o',markeredgecolor='orange', markerfacecolor='white',label=r'$f(y)= a y (b-y)$')
    #ax2.plot(country['Deltay'],country['preval'],'g.',label=r'$\partial_t y= a y (1-b-z-y)$')
    ax2.plot(country['incid']/country['popSize'],country['preval'],'.',color='blue',label=r'$(Y_{new,k} -     Y_{new,k-1}, y)$')
    #ax3.plot(country['dy'],country['preval'],'o',markeredgecolor='blue', markerfacecolor='white',label=r'$(\Delta Y_k , y)$')
    ax2.plot([0,0],[0,country['preval'].max()],'k:')
    ax2.plot([0,0],[0,country['preval'].max()],'k:')
    ax3.plot(days,country['alpha'],'.',color='blue',label=r'$\alpha$')
    tax3.plot(days,country['beta'],'.',color='orange',label=r'$\beta$')
    #ax1.set_ylim(0,1.1* country['preval'].max())
    #ax2.set_ylim(0,1.1* country['preval'].max())    
    #ax2.set_xlim(dy.min(),1.3*dy.max())
    ax1.set_xlabel(r'$t$ (days)');ax2.set_xlabel(r'$\partial_t y$'); ax3.set_xlabel(r'$t$ (days)');
    ax1.set_ylabel(r'$y_k$');ax2.set_ylabel(r'$y_k$')
    ax1.legend(loc='center left'); tax1.legend(loc='upper left')
    ax1.set_ylim(0,maxY); ax2.set_ylim(0,maxY); tax1.set_ylim(0,10); 
    ax3.set_ylim(0,10);     tax3.set_ylim(0,1)
    ax2.legend(loc='lower right'); 
    ax3.legend(loc='upper left'); tax3.legend(loc='center left'); 
    gr.ion(); gr.draw()
    fN= '../figures_COVID19_dataAnalysis/dam_COVID19_JHU_phasePlane_%s.png'%country['name']
    ff.suptitle('Phase plane, %s'%country['name'])
    ff.subplots_adjust(left=0.075,bottom=0.075,right=0.9,top=0.9,wspace=0.2,hspace=0.25)
    ff.savefig(fN)
    print('Saving file to %s'%fN)
    return ff


# Examples of phase plane dynamics

# In[499]:


Italy= getSingleCountryData(G, popSize= Pops_Millions, cou='Italy')
Italy= dataSIR_Country(country=Italy,a=230, b= 0.004)
fItaly= plotDataPhasePortraits(Italy)


# In[500]:


France= getSingleCountryData(G, popSize= Pops_Millions, cou='France')
France= dataSIR_Country(country=France,a=480, b= 0.002)
fFrance= plotDataPhasePortraits(France)


# In[509]:


US= getSingleCountryData(G, popSize= Pops_Millions, cou='US')
US= dataSIR_Country(US,a=100, b= 0.01)
fUS= plotDataPhasePortraits(US)


# In[515]:


Spain= getSingleCountryData(G, popSize= Pops_Millions, cou='Spain')
Spain= dataSIR_Country(Spain,a=370, b= 0.003)
fSpain= plotDataPhasePortraits(Spain)


# In[530]:


Belgium= getSingleCountryData(G, popSize= Pops_Millions, cou='Belgium')
Belgium= dataSIR_Country(Belgium,a=310, b= 0.0033)
fBelgium= plotDataPhasePortraits(Belgium)


# In[546]:


China= getSingleCountryData(G, popSize= Pops_Millions, cou='China')
China= dataSIR_Country(China,a=900, b= 0.001)
fChina= plotDataPhasePortraits(China)


# In[548]:


Mexico= getSingleCountryData(G, popSize= Pops_Millions, cou='Mexico')
Mexico= dataSIR_Country(Mexico,a=350, b= 0.0025)
fMexico= plotDataPhasePortraits(Mexico)


# In[ ]:





# In[ ]:




