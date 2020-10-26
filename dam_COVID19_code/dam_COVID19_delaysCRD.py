#!/usr/bin/env python
# coding: utf-8

# # Delays between case, recovery, and death reports during COVID-19 from different perspectives
# ## Marco Arieli Herrera-Valdez$^1$, Carlos Ignacio Herrera-Nolasco$^1$, 
# ## Alejandro Joel Herrera-McKiernan$^2$, Emilio Arieli Herrera-McKiernan$^3$, 
# ## Eugenia O'Reilly-Regueiro$^4$,
# 
# 
# #### $^1$ Departamento de Matemáticas, Facultad de Ciencias, Universidad Nacional Autónoma de México
# #### $^2$ Escuela Primaria República de Guatemala, Secretaría de Educación Pública, México
# #### $^3$ Escuela Secundaria Vicente Guerrero, Secretaría de Educación Pública, México
# #### $^4$ Instituto de Matemáticas, Facultad de Ciencias, Universidad Nacional Autónoma de México
# 
# Last modified: MAHV, 20200507
# 
# 
# 
# ### Data sources
# 
# The data was downloaded from the repository for the 2019 Novel Coronavirus Visual Dashboard operated by the Johns Hopkins University Center for Systems Science and Engineering (JHU CSSE) (https://github.com/CSSEGISandData/COVID-19). All calculations were performed using Python version 3.82 (https://www.python.org/) and the modules numpy (https://numpy.org/), matplotlib (https://matplotlib.org/), and pandas (https://pandas.pydata.org/). A JuPyTeR notebook with the analysis and calculations performed here can be found at (https://scab-unam.github.io/dam_COVID-19/).
# 

# In[1]:


from dateutil.parser import parse
import sys
sys.path.insert(1, '../')
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
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


def openCSV_DB(path,comp='zip',enc='latin-1'):
    data=pd.read_csv(path, compression=comp,encoding=enc)
    print('Data obtained from %s'%path)
    return data

sitio='https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
archCasos='time_series_covid19_confirmed_global.csv'
archFatal='time_series_covid19_deaths_global.csv'
archRecov='time_series_covid19_recovered_global.csv'

compr=None; codif='latin-1'
cases = openCSV_DB(sitio+archCasos,compr,codif)
recovs= openCSV_DB(enc=codif,path=sitio+archRecov,comp=compr)
deaths= openCSV_DB(path=sitio+archFatal,comp=compr)
#
locs_url='https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/ecdc/locations.csv'
locs= openCSV_DB(path= locs_url, comp=None)


# In[108]:


cases.head(10)


# In[110]:


cases.describe()


# In[3]:


cc='China';cases[cases['Country/Region'].isin([cc,])]


# In[4]:


locs.tail(5)


# In[109]:


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
days=np.arange(nDays)
print('Got data from %d days between %s and %s'%(nDays,dates[0],dates[-1]))
lastDay=parse(dates[-1])
print(lastDay.date())


# In[8]:


# -------------------
print("""Separating data by taking whole country cases into account""")
# -------------------
npCases = cases.to_numpy()
countries = np.unique(npCases[:,1])
nCountries = len(countries)
print('Considering data from {d} countries'.format(d=nCountries))
print(countries)


# Get a dictionary with data from each country

# In[9]:


cc='United Kingdom'
iC = cases['Country/Region'].isin([cc]); print(iC[iC==True])
cases[cases['Country/Region'].isin([cc,])]


# In[10]:


cc='Sweden'
iC = recovs['Country/Region'].isin([cc]); print(iC[iC==True])
recovs[recovs['Country/Region'].isin([cc,])]


# In[159]:


def diff(a):
    da = np.zeros(len(a))
    da[1:]= a[1:]-a[:-1]
    return da

def countryCases(countries, casesDF,recovDF,deathDF,locs):
    dd=dict()
    npCases = casesDF.to_numpy()
    npRecov = recovDF.to_numpy()
    npDeath = deathDF.to_numpy()
    countries = np.unique(npCases[:,1])
    nCountries = len(countries)
    dd['dates']= casesDF.columns[4:]
    for cc in countries:
        #print(cc)
        dd[cc] = dict()
        ii = locs['location'].isin([cc,])
        ps= locs[ii]['population'].to_numpy()
        dd[cc]['popSize'] =ps
        iC = cases['Country/Region'].isin([cc])
        dd[cc]['cCases']= cases[iC].iloc[:,4:].to_numpy().sum(0)
        iR = recovs['Country/Region'].isin([cc])
        dd[cc]['cRecov']= recovs[iR].iloc[:,4:].to_numpy().sum(0)
        iD = deaths['Country/Region'].isin([cc])
        dd[cc]['cDeath']= deaths[iD].iloc[:,4:].to_numpy().sum(0)
        nPts =len(dd[cc]['cCases'])
        #
        dd[cc]['pCases']= diff(dd[cc]['cCases'])
        # dd[cc]['pCases']= dd[cc]['cCases'] - dd[cc]['cRecov'] - dd[cc]['cDeath']
        dd[cc]['newRecov']= diff(dd[cc]['cRecov'])
        dd[cc]['newDeath']= diff(dd[cc]['cDeath'])
        dd[cc]['newRemovals'] = dd[cc]['newRecov']+dd[cc]['newDeath']
        dd[cc]['newCases'] = diff(dd[cc]['pCases']) + dd[cc]['newRemovals']
        dd[cc]['cCaseRemovalRatio'] = dd[cc]['newRemovals'] / dd[cc]['pCases']
        dd[cc]['recovRemovalRatio'] = dd[cc]['newRecov'] / dd[cc]['newRemovals']
        dd[cc]['deathRemovalRatio'] = dd[cc]['newDeath'] / dd[cc]['newRemovals']
        #dd[cc]['rateDeath'] = np.zeros(nPts); dd[cc]['rateRecov'] = np.zeros(nPts)
        #dd[cc]['timeDeath'] = np.zeros(nPts); dd[cc]['timeRecov'] = np.zeros(nPts)
        #ii = np.where(dd[cc]['pCases']>0)[0]
        #dd[cc]['rateDeath'][ii] = dd[cc]['newDeath'][ii]/dd[cc]['pCases'][ii]
        #dd[cc]['rateRecov'][ii] = dd[cc]['newRecov'][ii]/dd[cc]['pCases'][ii]
        #dd[cc]['timeDeath'][ii] = dd[cc]['pCases'][ii]/dd[cc]['newDeath'][ii]
        #dd[cc]['timeRecov'][ii] = dd[cc]['pCases'][ii]/dd[cc]['newRecov'][ii]
        jj = np.where(dd[cc]['newRecov']+dd[cc]['newDeath'] >0)[0]   
        dd[cc]['R(t)'] = np.zeros(nPts)
        dd[cc]['R(t)'][ii] = dd[cc]['newCases'][ii]/ (dd[cc]['newRecov'][ii]+dd[cc]['newDeath'][ii])
    return dd


# In[160]:


data= countryCases(countries, cases, recovs, deaths,locs)


# In[157]:


print('Keys for the global data dictionary: %s'%data['China'].keys())
cc = 'Germany'
print('Calculations for %s'%cc)
print('cRecov',data[cc]['cRecov'])
print('Recov times',data[cc]['timeRecov'])
print('Death times',data[cc]['timeDeath'])


# ### Delays between first case and first recovery or first death reports

# In[143]:


def findFirst(condition,defaultValue=None):
    i = np.where(condition)
    #print('len i',len(i[0]))
    if len(i[0])<1:
        j=defaultValue
    else: 
        j=i[0].min()
    return j

def firstCRDs(data):
    i1stC=list(); i1stD=list(); i1stR=list() 
    dfV=-1
    for nn in range(nCountries):
        cName=countries[nn]
        i1stC.append(findFirst(data[cName]['cCases']>0,defaultValue=dfV))
        i1stR.append(findFirst(data[cName]['cRecov']>0,defaultValue=dfV))
        i1stD.append(findFirst(data[cName]['cDeath']>0,defaultValue=dfV))
    return np.array(i1stC),np.array(i1stR),np.array(i1stD)

def checkedDelay(A,B):
    nPts = np.minimum(len(A),len(B))
    d=list()
    for n in range(nPts):
        if type(A[n])==type(B[n]):
            d.append(A[n]-B[n])
    return d


# #### Sorting the first case, first recovery, and first death reports for quantification

# In[144]:


# Indices of the first cases, recoveries, and deaths (these correspond to dates from d0), for each country
i1stC,i1stR,i1stD=firstCRDs(data)
# Sorting in ascending order for the occurrence of the first case reports. These indices correspond to countries, sorted by their date of first case report.
isort_i1stC = i1stC.argsort() 
# Sorting of the first report date by country using the sorted indices for the first cases 
si1stC= i1stC[isort_i1stC] 
si1stR= i1stR[isort_i1stC]
si1stD= i1stD[isort_i1stC]
# Difference between the indices for first recoveries and deaths relative to the first case reports (these correspond to days between first case report and first recovery or first death reports in each country)
delaysDC= i1stD-i1stC
# Delays by country sorted with respect to the first case reports
#sDelaysRC = delaysRC[isort_i1stC]
sDelaysDC = delaysDC[isort_i1stC]
# Countries listed 
s1stC_countries = countries[isort_i1stC] 


# ## Delays relative to the first case report in different countries
# 
# To get an idea of the initial the dynamics of the pandemic, we calculated the delays to the first case reports and compare to the first deaths reported by country. 

# #### Table to see countries by their first case report dates relative to $d_0$ and the delays to the first recovery or first death report

# In[145]:


# As many as countries
print('Country & First case & first recovery & first death \\\\')
print('\hline ')
fn= 'COVID19_worldFirstReports.tex'
file1 = open(fn,"w+")
for mm in range(nCountries):
    if ((si1stR[mm]>=0) & (si1stD[mm]>=0)):
        sss = [s1stC_countries[mm], dates[si1stC[mm]], si1stC[mm], dates[si1stR[mm]], si1stR[mm]-si1stC[mm], dates[si1stD[mm]], si1stD[mm]-si1stC[mm]]
        str0='{s[0]} & {s[1]} ({s[2]})  & {s[3]} ({s[4]})  & {s[5]} ({s[6]}) \\\\'.format(s=sss)
    elif ((si1stR[mm]<0)&(si1stD[mm]>=0)):
        sss = [s1stC_countries[mm], dates[si1stC[mm]], si1stC[mm], dates[si1stD[mm]], si1stD[mm]-si1stC[mm]]
        str0='{s[0]} & {s[1]} ({s[2]})  & -- (--)  & {s[3]} ({s[4]}) \\\\'.format(s=sss)
    elif (si1stD[mm]<0&(si1stR[mm]>=0)):
        sss = [s1stC_countries[mm], dates[si1stC[mm]], si1stC[mm], dates[si1stR[mm]], si1stR[mm]-si1stC[mm]]
        str0='{s[0]} & {s[1]} ({s[2]})  & {s[3]} ({s[4]})  & -- (--) \\\\'.format(s=sss)
    else:
        sss = [s1stC_countries[mm], dates[si1stC[mm]], si1stC[mm]]
        str0='{s[0]} & {s[1]} ({s[2]})  & -- (--)  & -- (--) \\\\'.format(s=sss)
    file1.write(str0); print(str0)
    
    if (mm<nCountries-1):
        if (si1stC[mm]<si1stC[mm+1]):
            file1.write('\hline'); print('\hline ')
file1.close()


# In[146]:


print('Days of first cases',i1stC)
print('Indices indices of the first cases sorted by date\n',isort_i1stC)


# In[147]:


print('Days of first death reports (not sorted)\n',i1stD)
print('Latest first recovery reported on day %d'%i1stR.max())
print('Latest first death reported on day %d'%i1stD.max())


# In[148]:


print('Delays between first case and first death reports in each country (not sorted)\n',delaysDC)


# #### Counting first case reports by day

# In[149]:


nDailyReports_C = np.zeros(nDays)
nDailyReports_R = np.zeros(nDays)
nDailyReports_D = np.zeros(nDays)
mm= np.maximum(nDays,si1stR.max())
#print(mm)
for nn in range(mm):
    nDailyReports_C[nn] = (si1stC==nn).sum()
    nDailyReports_R[nn] = (si1stR==nn).sum()
    nDailyReports_D[nn] = (si1stD==nn).sum()

print('# of first case reports as a function of the day:\n',nDailyReports_C)
print('# of first recovery reports as a function of the day:\n',nDailyReports_R)
print('# of first death reports as a function of the day:\n',nDailyReports_D)

cDailyReports_C = nDailyReports_C.cumsum()
cDailyReports_R = nDailyReports_R.cumsum()
cDailyReports_D = nDailyReports_D.cumsum()


# In[150]:


#cFirstCases, binsC = np.histogram(dFC,np.arange(0,(iFD-iFC).max()))
days = np.arange(nDays)

def plotCaseArrivals():
    shift=0.3; W=.7
    fRepArrivals= gr.figure(figsize=(11,11)); 
    gr.ioff(); rows=3; cols=1;
    ticks = np.arange(0,nDays,7)
    ax=list(); tAx=list()
    for m in range(cols*rows):
        ax.append(fRepArrivals.add_subplot(rows,cols,m+1))
        tAx.append(ax[m].twinx())
    tAx[0].plot(cDailyReports_C,'b',lw=2, label='Cumulative # countries reporting fist cases')
    tAx[1].plot(cDailyReports_C,'b',lw=2, label='Cumulative # countries reporting first cases')
    tAx[1].plot(cDailyReports_R,'orange',lw=2,label='Cumulative # countries reporting first deaths')
    tAx[2].plot(cDailyReports_C,'b',lw=2, label='Cumulative # countries reporting first cases')
    tAx[2].plot(cDailyReports_R,'orange',lw=2, label='Cumulative # countries reporting first recoveries')
    tAx[2].plot(cDailyReports_D,'k',lw=2,label='Cumulative # countries reporting first deaths')
    ax[0].bar(days,nDailyReports_C,color='blue',align='center',width=W,label='First case report')
    ax[1].bar(days,nDailyReports_R,color='orange',align='center',width=W,label='First recovery report')
    ax[2].bar(days,nDailyReports_D,color='black',align='center',width=W,label='First death report')
    for m in range(cols*rows):
        ax[m].set_xticks(ticks)
        tAx[m].set_yticks(np.arange(0,nCountries+1,10))
        ax[m].set_yticks(np.arange(0,nDailyReports_C.max()+1,2))
        ax[m].set_ylabel('# of countries',rotation=90)
        tAx[m].set_ylabel('# Cumulative # of countries',rotation=-90,labelpad=15)
        ax[m].legend(loc='upper left')
        tAx[m].legend(loc='center right')
    ax[2].set_xticklabels(dates[ticks],{'fontsize':8})
    for label in ax[2].get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('center')
        label.set_fontsize(8)
        label.set_horizontaloffset=-15
    for mm in range(rows*cols):
        ax[mm].set_xlabel('Days from $d_0$')
    fRepArrivals.subplots_adjust(left=0.075,bottom=0.1,right=0.9,top=0.95,wspace=0.1,hspace=0.25)
    fRepArrivals.suptitle('Delays in first case report, first death, and first recovery, all countries')
    gr.ion(); gr.draw(); gr.show()
    fRepArrivals.savefig('../figures_COVID19_dataAnalysis/dam_COVID19_JHU_reportArrivals_AllCountries_%s.png'%lastDay.date())
    return fRepArrivals


# In[151]:


fRepArrivals= plotCaseArrivals()


# In[152]:


# As many as countries
print('Country & First case & first recovery & first death \\\\')
print('\hline ')
iKeptCountries=list()
sDelay_1stRC=list()
sDelay_1stDC=list()
for mm in range(nCountries):
    if si1stD[mm]>=0:
        iKeptCountries.append(mm)
        sDelay_1stRC.append(si1stR[mm]-si1stC[mm])
        sDelay_1stDC.append(si1stD[mm]-si1stC[mm])
print(np.array(sDelay_1stRC))
print(sDelay_1stDC)


# In[210]:


#cFirstCases, binsC = np.histogram(dFC,np.arange(0,(iFD-iFC).max()))
days = np.arange(nDays)

def cumCaseArrivals(figS=(11,13)):
    shift=0.3; W=.7
    fRepArrivals= gr.figure(figsize=figS); 
    gr.ioff(); rows=2; cols=1;
    ticks = np.arange(0,nDays,7)
    ax=list(); tAx=list()
    d1=12; d2=19; d3=16; d4=20; L1=50; L2=140; 
    for m in range(cols*rows):
        ax.append(fRepArrivals.add_subplot(rows,cols,m+1))
    a = np.where(cDailyReports_C>=L1)[0][0]
    b = np.where(cDailyReports_C>=L2)[0][0]
    c = np.where(cDailyReports_R<L1)[0][-1]
    d = np.where(cDailyReports_R<L2)[0][-1]
    e = np.where(cDailyReports_D<L1)[0][-1]
    f = np.where(cDailyReports_D>=L2)[0][0]
    ax[0].plot([days[a],days[c]],[L1,L1],'k-',lw=2, alpha=0.35)
    ax[0].plot([days[b],days[d]],[L2,L2],'k-',lw=2, alpha=0.35)
    ax[0].plot(days,cDailyReports_R,'orange',lw=2,label='Cumulative # countries reporting deaths')
    ax[0].plot(days,cDailyReports_C,'b',lw=2, label='Cumulative # countries reporting cases')
    ax[0].plot(days[d1:],cDailyReports_C[:-d1],'b:',lw=1, alpha=0.5, label=r'$C(t+%d)$'%d1)
    ax[0].plot(days[d2:],cDailyReports_C[:-d2],'b--',lw=1, alpha=0.5, label=r'$C(t+%d)$'%d2)
    ax[1].plot([days[a],days[e]],[L1,L1],'k-',lw=2, alpha=0.35)
    ax[1].plot([days[b],days[f]],[L2,L2],'k-',lw=2, alpha=0.35)
    ax[1].plot(days,cDailyReports_C,'b',lw=2, label='Cumulative # countries reporting cases')
    ax[1].plot(days[d3:],cDailyReports_C[:-d3],'b:',lw=1, alpha=0.5, label=r'$C(t+%d)$'%d3)
    ax[1].plot(days[d4:],cDailyReports_C[:-d4],'b--',lw=1, alpha=0.5, label=r'$C(t+%d)$'%d4)
    ax[1].plot(days,cDailyReports_D,'k',lw=2,label='Cumulative # countries reporting deaths')
    for m in range(cols*rows):
        ax[m].set_xticks(ticks)
        ax[m].set_yticks(np.arange(0,nCountries+1,10))
        ax[m].set_ylabel('# Countries',rotation=90,labelpad=15)
        ax[m].legend(loc='lower right',fontsize=13)
    ax[1].set_xticklabels(dates[ticks],{'fontsize':8,'rotation':45})
    for label in ax[1].get_xticklabels():
        label.set_horizontalalignment('center')
        label.set_horizontaloffset=-15
    for mm in range(rows*cols):
        ax[mm].set_xlabel('Days from $d_0$')
    fRepArrivals.subplots_adjust(left=0.075,bottom=0.1,right=0.9,top=0.95,wspace=0.1,hspace=0.25)
    fRepArrivals.suptitle('Delays in first case report, first death, and first recovery, all countries')
    gr.ion(); gr.draw(); gr.show()
    str1='../figures_COVID19_dataAnalysis/dam_COVID19_JHU_firstDelays_AllCountries_%s.png'%lastDay.date()
    fRepArrivals.savefig(str1)
    return fRepArrivals


# In[211]:


fRepArrivals= cumCaseArrivals(figS=(11,13))


# ## Prevalence and disease dynamics at the country level

# ### Simple model of epidemiological dynamics based on macroscopic data
# 
# $C(t)$ total confirmed cases, $R(t)$ total recoveries (not always available), $D(t)$ total deaths (not always available)
# 
# $\delta$ is the time step in days. 
# Absolute prevalence (active cases without normalization by population size)
# $$c(t + \delta) = C(t+\delta) - C(t)$$ 
# If $\delta =1$, then $c$ represents the daily prevalence (without normalization by population size)
# 
# 
# 
# 
# 
# 

# #### Recoveries and deaths 
# The new recoveries and deaths can be calculated by substraction
# \begin{eqnarray}
# r(t+\delta) &=& R(t+\delta) - R(t) \\
# d(t+\delta) &=& D(t+\delta) - D(t) 
# \end{eqnarray}
# 
# Incidence (new cases) ($x$)
# $$ x(t) = \frac{c(t+\delta) - c(t)}{\delta} + r(t) + d(t) $$
# 
# 
# The proportion of new recoveries and deaths from the existing cases is
# $$ \rho = \frac{r(t) + d(t)}{c(t)}$$
# The event in which an active case either recovers or dies (gets removed from the active cases) can be assumed to be independent of other cases. 
# Then the number of new removals from the case population can be assumed to be a binomial random variable with parameters $c(t)$ and $\rho(t)$, and $\rho(t)\cdot c(t)$ can be assumed to be the expected number of recoveries+deaths. 
# 
# The proportion of deaths from the total of deaths and recoveries at time t, $p_d$ can be regarded as the probability that a removal results in death. Similarly, the proportion of recoveries within the removals, $p_r$ can be regarded as the probability that a removal is a recovery.
# 
# Explicitly
# $$ p_r = \frac{r(t)}{r(t) + d(t)}, \quad p_d = \frac{d(t)}{r(t) + d(t)}.$$
# 
# It is worth noticing that the deaths and recoveries occur with a delay with respect to the first date of positivity. 
# 
# 

# In[307]:


def plotEpiDynamics(cc = 'Germany'):
    nDays = len(days)
    cCases_norm0=data[cc]['cCases']/data[cc]['cCases'].max()
    cRecov_norm0=data[cc]['cRecov']/data[cc]['cCases'].max()
    cDeath_norm0=data[cc]['cDeath']/data[cc]['cCases'].max()
    pCases =  diff(data[cc]['cCases'])
    # print(cCases_norm0)
    f=gr.figure(figsize=(17,13)); rows=3; cols=1; gr.ioff(); 
    ax=list(); 
    ax.append(f.add_subplot(rows,cols,1))
    ax[0].plot(days,cCases_norm0,lw=2,label=r'$C(t)/\max(C(t))$ Normalized total confirmed cases')
    ax[0].plot(days,cRecov_norm0,label=r'$R(t)/\max(C(t))$, Normalized recoveries')
    ax[0].plot(days,cDeath_norm0,label=r'$D(t)/\max(C(t))$, Normalized total deaths')
    ax.append(f.add_subplot(rows,cols,2))
    ax[1].plot(days,data[cc]['pCases'],lw=3,label=r'(absolute) prevalence')
    ax[1].plot(days,data[cc]['newRecov'],label=r'recoveries')
    ax[1].plot(days,data[cc]['newDeath'], label=r'deaths')
    ax[1].set_ylim(0,data[cc]['pCases'].max())
    ax.append(f.add_subplot(rows,cols,3))
    ax[2].plot(days,data[cc]['cCaseRemovalRatio'],'.',label=r'$\rho(t)$')
    ax[2].plot(days,data[cc]['recovRemovalRatio'],'.',label=r'$p_r(t)$')
    ax[2].plot(days,data[cc]['deathRemovalRatio'],'.',label=r'$p_d(t)$')
    ax[2].set_ylim(0,1)
    #
    ticks = np.arange(0,nDays,7)
    for nn in range(len(ax)):
        ax[nn].legend(fontsize=15)
        ax[nn].set_xticks(ticks)
        ax[nn].set_xticklabels(dates[ticks],{'fontsize':10,'rotation':45})
        for label in ax[nn].get_xticklabels():
            label.set_horizontalalignment('center')
            label.set_horizontaloffset=-15
    gr.ion(); gr.draw(); 
    str1='../figures_COVID19_dataAnalysis/dam_COVID19_JHU_prevalenceDynamics_%s_%s.png'%(cc,lastDay.date())
    fRepArrivals.savefig(str1)
    #meanCFR = data[cc]['deathRemovalRatio'].mean()*100
    #print('Mean case fatality ratio = %g per cent'%(meanCFR))
    return ff,ax


# In[308]:


ccc='Japan'
plotEpiDynamics(cc = ccc)


# ## Case fatality ratios
# 
# The case fatality ratio calculation for an ongoing epidemic should be $100\cdot p_d$ [(Estimating mortality during a pandemic, WHO, 2020)](https://www.who.int/news-room/commentaries/detail/estimating-mortality-from-covid-19). For that, it is necessary to assume that 
# 
# 1. The likelihood of detecting cases and deaths is consistent over the course of the outbreak.
# 
# 2. All detected cases have resolved (that is, reported cases have either recovered or died).
# 
# Therefore, the data to take into consideration must be such that the values of the distribution of $p_d$ are relatively stable. 

# In[306]:


upperLim=0.7
cfr0 = data[ccc]['cDeath']/(data[ccc]['cDeath']+data[ccc]['cRecov'])
ii = (cfr0<upperLim) & (cfr0>0)
cfr = cfr0[ii]
#
fCFR=gr.figure(figsize=(17,5)); gr.ioff()
ax=list(); rows = 1; cols=2;
ax.append(fCFR.add_subplot(rows,cols,1))
ax.append(fCFR.add_subplot(rows,cols,2))
ax[0].plot(days[ii],cfr,'.')
ax[1].hist(cfr,np.arange(0,1,0.01),orientation='horizontal')
ax[1].set_xlabel(ccc+ ' CFRs' );
ax[1].set_yticks(np.arange(0,0.5,0.1))
for i in range(rows*cols):
    ax[i].set_ylim(0,np.maximum(upperLim,cfr.max()))
ticks = np.arange(0,nDays,14)
ax[0].legend(fontsize=15)
ax[0].set_xticks(ticks)
ax[0].set_xticklabels(dates[ticks],{'fontsize':10,'rotation':45})
for label in ax[0].get_xticklabels():
    label.set_horizontalalignment('center')
    label.set_horizontaloffset=-15

gr.ion(); gr.draw()


# ## Deterministic epidemic dynamics by country from the simple model 
# 
# For a start, assume that the population of each country is divided into four subsets representing the susceptible, infected, no-longer infected (non-infectious) and no-longer susceptible, and dead due to infection, with densities respectively written as 
# $u$, $v$, $r$, and $d$. Assume that $1=u+v+w+x$, and that the dynamics between those groups follow an SIR-like evolution rules given by
# \begin{eqnarray}
# \partial_t u &=& -\lambda u ,\\
# \partial_t v &=& \lambda u - v \left( \frac{p_r}{\tau_r} +\frac{1-p_r}{\tau_d} \right) ,\\
# \partial_t r &=& v\frac{p_r}{\tau_r}, \\
# \partial_t d &=& v\frac{1-p_r}{\tau_d},
# \end{eqnarray}
# with $\tau_r$ and $\tau_d$ representing the average times to clear the infection and, alternatively, the expected time to death due to disease. Although these two times could differ depending on a number of factors including the health care capacities in different countries, they may not be very different between different countries.  It is not unreasonable to assume that is the case to construct an initial model that describes the macroscopic dynamics of the COVID-19 epidemics at the whole country level.   Note the desinfection rate and death rate are $p_r/\tau_r$ and $(1-p_r)/\tau_x$ respectively, with $p_r$ representing the proportion of recoveries among removals. 
# 
# 
# 

# In[269]:


def plotDynamicsPP(cc,figS=(15,7)):
    ff= gr.figure(figsize=figS)
    rows=1; cols=2; gr.ioff()
    ff.suptitle(cc)
    ax=list()
    for i in range(rows*cols):
        ax.append(ff.add_subplot(rows,cols,i+1))
              
    ax[0].plot(data[cc]['pCases'],'k',alpha=0.35, ms=2, label=r'Cases')
    ax[0].plot(data[cc]['newCases'],'.',ms=5,label=r'New cases')
    ax[0].plot(data[cc]['newRecov'],'.',label=r'recoveries')
    ax[0].plot(data[cc]['newDeath'],'.',label=r'deaths')
    ax[0].set_ylim(-1,data[cc]['newCases'].max())
    ax[0].set_ylabel('# People')
    ax[1].plot(data[cc]['pCases'], data[cc]['newCases'],'.',ms=5,label=r'New cases')
    ax[1].plot(data[cc]['pCases'], data[cc]['newRecov'],'.',label=r'recoveries')
    ax[1].plot(data[cc]['pCases'], data[cc]['newDeath'],'.',label=r'deaths')
    ax[1].set_ylim(-1,data[cc]['newCases'].max())
    ax[1].set_xlabel('Cases')
    ax[1].set_ylabel('New cases')
    for i in range(rows*cols):
        ax[i].legend(fontsize=15)
    return ff,ax


# In[270]:


ff,ax= plotDynamicsPP(cc = 'Belgium',figS=(17,7))


# The _force of infection_ in the model is given by $\lambda = \epsilon \alpha v$ where $\alpha$ is the rate of infection given an infectious contact, and $\epsilon$ is an exposure factor that scales the number of contacts between susceptibles and infected individuals. 
# As a consequence, it is possible to analyze the macroscopic dynamics with two equations in mind,
# \begin{eqnarray}
# \partial_t v &=& v \left[\alpha \epsilon \left(1-w-x-v\right)- \left(\frac{p}{\tau_w} +\frac{1-p}{\tau_x}\right) \right]  ,\\
# \partial_t (w+x) &=& v \left( \frac{p}{\tau_w} +\frac{1-p}{\tau_x} \right).
# \end{eqnarray}
# From there, the condition for epidemic dynamics is that 
# $$
# R(t)=\frac{\alpha \epsilon}{\left( \frac{p}{\tau_w} +\frac{1-p}{\tau_x} \right)} (1-w(t)-x(t)-v(t))> 1.  
# $$
# As a consequence, 
# $$R_o \approx \frac{\alpha \epsilon}{\frac{p}{\tau_w} +\frac{1-p}{\tau_x}}$$
# 
# 
# 

# Data from different countries suggest that on average, a person spends between 5 and 7 days in a stage that can be assumed to be primary viremia, only allowing virus replication, but infectious since the first hours after contagion. Then some people develop the disease, 50% (IQR) of which do so between 5 and 11 days. Importantly, people can be regarded as infectious from the day they are in first contact with SARS-CoV-2 \citep{}, and remain infectious until after clinical recovery.  From there, people who have been tested positive for SARS-CoV-2 with RT-PCR testing have been reported to remain infectious between 15 and 21 more days, depending on whether they recover or whether they die. This means that the __waiting time for people that recover from having COVID-19 can remain infectious approximately between 20 and 22 days__. Similarly, __people that die due to COVID-19 may remain infected between 26 and 28 days__. At the time this report is written, _there have been some reports claiming that people who have been tested positive for COVID-19 may still be positive for SARS-CoV-2 tests up to 160 days after the first time they tested positive_. However, that evidence is not conclusive yet.
# 
# According to data published before or on May 25, 2020, _the waiting times for removal from the infectious group can be assumed to be aproximately 20 days for $\tau_w$ and 25 days for $\tau_x$_. 
# The __initial estimates for $R_o$ in different countries fluctuated between 2 and 4__. A crude estimate for the product of the exposure factor and the infection rate given an infectious contact, can be obtained by using the estimates for $R_o$. First, assume that at the time of the outbreak, $1-v-w-x \approx 1$. Therefore, $$\frac{p}{\tau_w} +\frac{1-p}{\tau_x} =\frac{p}{20} +\frac{1-p}{25} = \frac{25(1-p) + 20p}{500}$$ and
# $$
# \alpha \epsilon = R_o \left( \frac{p_r}{\tau_r} +\frac{1-p_r}{\tau_d} \right) = \left(\frac{25(1-p) + 20p}{500}\right)R_o  
# $$
# A static, conservative estimate can be obtained assuming that the case fatality ratio is similar to $p_d = 1-p_r$, which is between 0.005 and 0.0158 [(Health Metrix)](www.healthdata.org). 

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# The data available for most countries is separated in three cumulative time series, respectively, the total cases ($C$), recoveries ($R$), and deaths ($D$), sampled every day. After normalization by the population size $T$, recoveries and deaths can be thought of as similar in behavior to $w+x$. The prevalence of the disease cannot be known. However, if testing is done in a way that the sampling is representative of the overall dynamics, an approximation with a similar behaviour to the prevalence can be obtained from the the data of confirmed cases is then
# $$
# \tilde{P}(t)=\frac{C(t) - \left( R(t) + D(t) \right)}{T}
# $$
# which means that the dynamics of $v(t)$ can be fit to the dynamics of $\tilde{P}$, in an attempt to study the qualitative features of the epidemic. The possible undersampling due to incomplete testing can then be accounted non explicitly by $\epsilon$. 
# Then, from the data, it is also possible to estimate the incidence and the rates. The change in the deaths per unit time could be assumed to be $$\partial_t x(t+h) \approx D(t+h) - D(t).$$ Similarly for the change in the recoveries $$\partial_t w(t+h) \approx R(t+h) -R(t).$$ 
# The death rate can be estimated by using the equation
# $$
# \frac{1-p}{\tau_x} \approx \frac{\partial_t x}{v} 
# $$
# The desinfection rate can be obtained in a similar way
# $$
# \frac{p}{\tau_w} \approx \frac{\partial_t w}{v} 
# $$
# 
# 

# In[139]:


# create a new plot with a title and axis labels
from bokeh.plotting import figure, show
from bokeh.io import output_notebook
output_notebook()

p = figure(title="simple line example", x_axis_label='x', y_axis_label='y')
cc = 'Canada'
# add a line renderer with legend and line thickness
p.line(days,data[cc]['cCases'], legend_label="Cases ", line_width=2,line_color='blue')
p.line(days,data[cc]['cRecov'], legend_label="Recoveries", line_width=2,line_color='green')
p.line(days,data[cc]['cDeath'], legend_label="Deaths", line_width=2,line_color='red')
p.line(days,data[cc]['pCases'], legend_label="prevalence", line_width=2,line_color='black')
# show the results
show(p)


# In[ ]:




