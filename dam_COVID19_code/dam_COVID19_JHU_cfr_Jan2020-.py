#!/usr/bin/env python
# coding: utf-8

# # Initial mortality and recovery in reported cases of COVID-19 a few weeks into the pandemic
# 
# ## Carlos Ignacio Herrera-Nolasco$^1$, Alejandro Joel Herrera-McKiernan$^2$, 
# ## Emilio Arieli Herrera-McKiernan$^3$, Eugenia O'Reilly-Regueiro$^4$,
# ## Marco Arieli Herrera-Valdez$^1$
# 
# #### $^1$ Departamento de Matemáticas, Facultad de Ciencias, Universidad Nacional Autónoma de México
# #### $^2$ Escuela Primaria República de Guatemala, Secretaría de Educación Pública, México
# #### $^3$ Escuela Secundaria Vicente Guerrero, Secretaría de Educación Pública, México
# #### $^4$ Instituto de Matemáticas, Facultad de Ciencias, Universidad Nacional Autónoma de México
# 
# Last modified: MAHV, 20200424
# 
# ### Abstract
# The time series of cases, recoveries, and deaths during a pandemic are usually contaminated by different problems including under report,  and delays that in many cases can be justified by late reactions by government agencies. 
# However, it should be possible to find general characteristics, and also important differences between the dynamics observed in different countries. Both could be important macroscopic markers to understand the global dynamics as a function of the local (country) dynamics of the disease caused by the pandemic pathogen. 
# In an attempt to extract general information from the global dynamics of COVID-19 cases, we analyzed the time series of cases, recoveries, and deaths reported up until April 24, 2020, and focused on delays between different types of reports, with special attention to the relationships between cases and fatalities.
# The first countries to report cases also reported fatalities within short or no delay. However, as more countries reported cases, consistent delay patterns showing a relatively constant delay between cases and recoveres, and between cases and fatalities within cases.  We also found a consistent contribution to of the case fatalities by age in countries like China, South Korea, and Italy, which have different health care systems, different mitigation capacities, and different population structures. We also examine the long term behaviour of the whole country CFR calculated from the time series of cases and deaths, both cumulative and from incidence, taking possible deaths-cases delays  between 0 and 30 days. The global data shows a strong indication of an asymptotic convergence toward a small interval of values in the delayed quotients, both cumulative and from incidence, in most countries that have at least 6 weeks reporting cases. For single countries, the values of the CFRs with delays between 0 and 30 days  for the last few days show an increasing trend as a function of the delay, but still between 0.05 and 0.3 for most countries with enough reports. Although the converging trend will not be clear until more data is available, the trend  shows that CFRs calculated by simply dividing deaths by cases without taking delays into account from macroscopic time series are not as incorrect as some reports indicate, at least with respect to delay, and without taking other biases like under-reporting into account.  
# 
# 
# ### Methods and data sources
# 
# The data was downloaded from the repository for the 2019 Novel Coronavirus Visual Dashboard operated by the Johns Hopkins University Center for Systems Science and Engineering (JHU CSSE) (https://github.com/CSSEGISandData/COVID-19). All calculations were performed using Python version 3.82 (https://www.python.org/) and the modules numpy (https://numpy.org/), matplotlib (https://matplotlib.org/), and pandas (https://pandas.pydata.org/). A JuPyTeR notebook with the analysis and calculations performed here can be found at (https://scab-unam.github.io/dam_COVID-19/).
# 
# 

# Module imports (all basic functions and other module specs can be found in [dam_COVID19_baseCode.py](dam_COVID19_baseCode.py)

# In[1]:


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
get_ipython().magic(u'matplotlib inline')


# In[2]:


cases, deathCases,recovCases = getCSSEGISandData(urlData=1)


# Let us first describe the data set. The data has been transposed so that the cases, deaths, and recovered the different regions of the world are stored in columns, each with the data from different dates starting from January 22, 2020. 

# In[3]:


cases.head(6)


# In[4]:


deathCases.tail(4)


# This means there are 4 columns before the time series begins, and the data contains 248 regions from the world. The data that can be used for calculations involving the cases in the epidemic can be found from the fifth row on.  Since the dates in the original data have only two numbers for the year, let us create a date range with the dates formatted with a four figures, and use the new dates later for the plots and illustrations.

# In[5]:


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
print('Got data from %d days between %s and %s'%(nDays,dates[0],dates[-1]))


# In[6]:


# -------------------
print("""""")
# -------------------
npCases = cases.to_numpy()
countries = np.unique(npCases[:,1])
nCountries = len(countries)
print('Considering data from {d} countries'.format(d=nCountries))


# Data curation observations.
# - Some countries report by province. This allows examination of the time series with a better spatial resolution.
# The countries where there are reports by province are China, United Kingdom, Australia, and Canada.
# - Arrays that include whole country cases, deaths, and recovered obtained
# and separated for analysis
# - Recovered are treated separately because reports do not include as many places as cases or deaths
# - Case fatality ratios are calculated to analyse the epidemia in different provinces and contrast with whole country estimates
# 

# In[7]:


# -------------------
# Gather same country data
# -------------------
gCases, countries_Cases= sortDataByCountry(cases, nHeaderCols)
gDeaths,countries_DeathCases = sortDataByCountry(deathCases, nHeaderCols)
gRecovCases,countries_RecovCases = sortDataByCountry(recovCases, nHeaderCols)

# -------------------
# Sum the counts from each country and construct a new array
# -------------------
# These arrays have the same size as the countries array (unique countries)
totCases=gatherDataByCountry(df=cases,nHeaderCols=4)
totDeathCases=gatherDataByCountry(df=deathCases,nHeaderCols=4)
totRecovCases=gatherDataByCountry(df=recovCases,nHeaderCols=4)


# ## Delays relative to the first case report in different countries
# 
# To get an idea of the initial the dynamics of the pandemic, we calculated the delays to the first case reports and compare to the first deaths reported by country. 

# In[8]:


# Describe the order of appearance of cases
dFC=findFirstCaseDates(totCases)
dFD=findFirstCaseDates(totDeathCases)
dFR=findFirstCaseDates(totRecovCases)
iArrival = dFC.argsort()
for n in range(37):
    print(countries[n],dFC[n],dFD[n],dFR[n],'rank: %dth'%iArrival[n])
 
d0SortedCountries = countries[iArrival]
sdFC = dFC[iArrival]
sdFD = dFD[iArrival]
sdFR = dFR[iArrival]


# In[9]:


# As many as countries
print('%First case, & to first death & to first recovery \\\\')
print('\hline &&&\\\\')
for mm in range(nCountries):
    #print('%s & %s & %s & %s\\\\'%(countries[jj[mm]], dates[siC[mm]], dates[siD[mm]], dates[siR[mm]]))
    #print('& %d & %d & %d \\\\ '%(siC[mm],siD[mm]-siC[mm],siR[mm]-siC[mm]))
    sss = [d0SortedCountries[mm], dates[sdFC[mm]], sdFC[mm], dates[sdFD[mm]], sdFD[mm]-sdFC[mm], dates[sdFR[mm]], sdFR[mm]-sdFC[mm]]
    print('{s[0]} & {s[1]} ({s[2]})  & {s[3]} ({s[4]})  & {s[5]} ({s[6]}) \\\\'.format(s=sss))
    if (mm<nCountries-1):
        if (sdFC[mm]<sdFC[mm+1]):
            print('\hline ')


# ## The time course of reports

# In[10]:


nCountriesFCByDate= list()
nCountriesFDByDate= list()
nCountriesFRByDate= list()
for m in range(nDays):
    # Find the indices of those who reported their first case on day m
    nCountriesFCByDate.append(len(np.where(dFC==m)[0])) 
    # Find the indices of those who reported their first death on day m
    nCountriesFDByDate.append(len(np.where(dFD==m)[0])) 
    nCountriesFRByDate.append(len(np.where(dFR==m)[0])) 
nCountriesFCByDate = np.array(nCountriesFCByDate)
nCountriesFDByDate = np.array(nCountriesFDByDate)
nCountriesFRByDate = np.array(nCountriesFRByDate)
cumC=nCountriesFCByDate.cumsum()
cumD=nCountriesFDByDate.cumsum()
cumR=nCountriesFRByDate.cumsum()


# In[11]:


print(nCountriesFCByDate,len(nCountriesFCByDate))
print(nCountriesFDByDate,len(nCountriesFDByDate))
print(nCountriesFRByDate,len(nCountriesFRByDate))


# The first countries that reported cases of COVID-19 on January 22, 2020 ($d_0$) were South Korea, China, Taiwan, US, Japan, and Thailand. Australia reported their first case 4 days later, together with Canada, on January 26, 2020. The first reports from Europe came 5 days later from Germany, on January 27. The first report in the Middle East came from the United Arab Emirates on January 29, 7 days after the first case reported. 
# Therefore, the COVID-19 epidemic had been reported from all continents but Africa one week after the first case was reported, and the first report from Africa arrived 23 days later from Egypt, on February 14, 2020.
# These data indicate that the COVID-19 epidemic was wide spread throughout the globe within 2 weeks of the first case, strongly suggesting that the dissemination of the SARS-CoV-2 virus among the human population occurred with a delay of approximately 2 weeks. 
# 
# These data reflect the delay in reports from different countries relative to $d_0$. Of note, 25 countries reported cases within 10 days from $d_0$, and 31 countries had reports within the first 30 days. 
# Of interest to the authors, Mexico reported their first case 37 days later, on February 28, 2020, together with Nigeria, New Zealand, Lithuania, Iceland, and Belarous. The last country that reported cases was Sao Tome and Principe, on April 6th, 2020, 75 days after the first report.
# 

# ### Delays between first case, first recovery, and first death reports

# In[12]:


#cFirstCases, binsC = np.histogram(dFC,np.arange(0,(iFD-iFC).max()))
days = np.arange(nDays)

def plotCaseArrivals():
    shift=0.3; W=.7
    fRepArrivals= gr.figure(figsize=(7,9)); 
    gr.ioff(); rows=3; cols=1;
    ticks = np.arange(0,nDays,7)
    ax=list(); tAx=list()
    for m in range(cols*rows):
        ax.append(fRepArrivals.add_subplot(rows,cols,m+1))
        tAx.append(ax[m].twinx())
    tAx[0].plot(cumC,'b',lw=2, label='Cumulative # countries reporting cases')
    tAx[1].plot(cumC,'b',lw=2, label='Cumulative # countries reporting cases')
    tAx[1].plot(cumD,'k',lw=2,label='Cumulative # countries reporting deaths')
    tAx[2].plot(cumC,'b',lw=2, label='Cumulative # countries reporting cases')
    tAx[2].plot(cumD,'k',lw=2,label='Cumulative # countries reporting deaths')
    tAx[2].plot(cumR,'orange',lw=2, label='Cumulative # countries reporting recoveries')
    ax[0].bar(days,nCountriesFCByDate,color='blue',align='center',width=W,label='First case report')
    ax[2].bar(days,nCountriesFDByDate,color='black',align='center',width=W,label='First death report')
    ax[1].bar(days,nCountriesFRByDate,color='orange',align='center',width=W,label='First recovery report')
    for m in range(cols*rows):
        ax[m].set_xticks(ticks)
        tAx[m].set_yticks(np.arange(0,nCountries+1,10))
        ax[m].set_yticks(np.arange(0,nCountriesFCByDate.max()+1,2))
        ax[m].set_ylabel('# Countries',rotation=90)
        tAx[m].set_ylabel('# Countries',rotation=-90,labelpad=15)
        ax[m].legend(loc='center left')
        tAx[m].legend(loc='upper left')
    ax[2].set_xticklabels(dates[ticks],{'fontsize':8})
    for label in ax[2].get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('center')
        label.set_fontsize(8)
        label.set_horizontaloffset=-15
    for mm in range(rows*cols):
        ax[mm].set_xlabel('Days from $d_0$')
    fRepArrivals.subplots_adjust(left=0.075,bottom=0.075,right=0.9,top=0.95,wspace=0.1,hspace=0.25)
    fRepArrivals.suptitle('Delays in first case report, first death, and first recovery, all countries')
    gr.ion(); gr.draw(); gr.show()
    fRepArrivals.savefig('./figures_COVID19_dataAnalysis/dam_COVID19_JHU_reportArrivals_AllCountries.png')
    return fRepArrivals

fRepArrivals= plotCaseArrivals()


# Approximately around 30 days after $d_0$, a relatively constant delay between first cases reported and the corresponding first recoveries and first deaths reported, lasting approximately 40 days. The biases in reports show up again after day 70 from $d_0$. In the beginning of the pandemic the biases are caused by short or no-delays between first cases and deaths. In the last portion of the period under consideration, the biases are caused by long delays in the first case, first recoveries, and first deaths.  
# The large bar in the recoveries from February 12, 2020, is caused by a large adjustment in the reports from China and the first recoveries reported in other countries. The large bar from the deaths at on April 23, 2020, is caused by late reporting from Africa. 

# In[13]:


n=7; print(r'%d countries reported within the first %d days from $d_0$'%(len(sdFC[sdFC<n+1]),n))
countries[iArrival[(sdFD-sdFC) <n+1]]


# The difference between the curve of cumulative number of countries reporting cases and the cumulative number of countries reporting deaths is very similar for a large interval of days. 
# The descriptive statistics of delay between death and case reports are as follows. 

# The delays between the first case reported and the first death reported within a single country have a wide distribution, between 0 and 76 days. 
# However, the bulk of the distribution can be found between 0 and 30 days, with a mode at 16 days. In consideration of the time course of the infection by SARS-CoV-2 reported in the literature \citep{}, the delay between death and contagion can be 21 or more days \citep{}. 
# 
# Assuming that the testing started at least one week before the date of the first case reported, the data suggest that the epidemic in those countries contributing to left portion of the histogram had started at least 3 weeks before their report. The mode at 16 days fits within the 5-7 days incubation periods reported plus the 2 weeks to develop severe symptoms. For those countries that reported their first death at the mode or after, these data suggest that the epidemic started in those countries at least one week before their first case reported.
# 
# The delays between the first recovery and the first case reported within a single country are distributed between 0 and 38 days. However, for most countries this delay is between 0 and 27 days approximately, with two clear modes at 12 and 16 days.  

# ### Simple estimation of the temporal bias induced by delays between first case, and first death or first recovery reports

# In[14]:


def findHDistance(t,x,y,m,M,nPts=100):
    if (len(t)!=len(x)) | (len(x)!=len(y) ):
        return 
    v = np.linspace(m,M,nPts)
    d=list()
    for n in range(len(v)):
        a = np.where(x>v[n])[0][0]
        b = np.where(y>v[n])[0][0]
        d.append(b-a)
    return np.array(d)

def findMode(x):
    xbins = np.unique(x)
    counts,bins =np.histogram(x,xbins)
    return x[counts.argmax()] 


# There is an interval where the delays are relatively constant. 

# In[15]:


v1 = 40; v2=150
delaysCumFDC=findHDistance(days,cumC,cumD,v1,v2,nPts=v2-v1)
delaysCumFRC=findHDistance(days,cumC,cumR,v1,v2,nPts=v2-v1)
print(delaysCumFDC)
print(delaysCumFRC)
mm = findMode(delaysCumFDC)
MM = findMode(delaysCumFRC)
print('The most common delay between the cumulative reports of cases and deaths is %d'%mm)
print('The most common delay between the cumulative reports of cases and recoveries is %d'%MM)


# In[16]:


# 
a = 50; b = 136
ia_cumC = np.where(cumC>a)[0].min()
ia_cumD = np.where(cumD>a)[0].min()-1
ib_cumC = np.where(cumC>b)[0].min()
ib_cumD = np.where(cumD>=b)[0].min()
ia_cumR = np.where(cumR>a)[0].min()-1
ib_cumR = np.where(cumR>=b)[0].min()
aa=ia_cumD-ia_cumC
bb=ib_cumD-ib_cumC
cc=ia_cumR-ia_cumC
dd=ib_cumR-ib_cumC
print(ib_cumC, ib_cumR, cc,dd)
print(ia_cumC,ib_cumC,ia_cumD,ib_cumD)
delaysFDC =sdFD-sdFC
delaysFRC =sdFR-sdFC
#
fRepDelay= gr.figure(figsize=(7,9)); gr.ioff()
rows=2; cols=2
ax=fRepDelay.add_subplot(rows,cols,1); 
bx=fRepDelay.add_subplot(rows,cols,2); 
cx=fRepDelay.add_subplot(rows,cols,3); 
dx=fRepDelay.add_subplot(rows,cols,4); 
days = np.arange(nDays)
cC = np.arange(nCountries)
ticks= np.arange(0,nDays,7)
bx.plot(cC, delaysFDC, '.', label='#case-death ')
bx.plot(cC, mm*np.ones(nCountries),label=r'%d days'%delaysFDC.mean())
dx.plot(cC, delaysFRC, '.', label='#case-death ')
dx.plot(cC, MM*np.ones(nCountries),label=r'%d days'%delaysFRC.mean())
ax.plot(days,cumC,label='Cases')
ax.plot(days,cumD, label='Deaths')
cx.plot(days,cumC,label='Cases')
cx.plot(days,cumR, label='Recoveries')
ax.plot(days,a*np.ones(nDays),':',color='orange',alpha=0.35)
ax.plot(days,b*np.ones(nDays),':',color='orange',alpha=0.35)
cx.plot(days,a*np.ones(nDays),':',color='orange',alpha=0.35)
cx.plot(days,b*np.ones(nDays),':',color='orange',alpha=0.35)
ax.plot( [days[ia_cumC],days[ia_cumC]], [cC[0],cC[-1]],':',color='orange',lw=1,alpha=0.35)
ax.plot( [days[ia_cumD],days[ia_cumD]], [cC[0],cC[-1]],':',color='orange',lw=1,alpha=0.35)
ax.plot( days[ib_cumC:ib_cumD], b*np.ones(bb),'k',lw=3,alpha=0.35)
ax.plot( days[ia_cumC:ia_cumD], a*np.ones(aa),'k',lw=3,alpha=0.35)
cx.plot( [days[ia_cumC],days[ia_cumC]], [cC[0],cC[-1]],':',color='orange',lw=1,alpha=0.35)
cx.plot( [days[ia_cumR],days[ia_cumR]], [cC[0],cC[-1]],':',color='orange',lw=1,alpha=0.35)
cx.plot( days[ib_cumC:ib_cumR], b*np.ones(dd),'k',lw=3,alpha=0.35)
cx.plot( days[ia_cumC:ia_cumR], a*np.ones(cc),'k',lw=3,alpha=0.35)
str1='Delay between first death and first case reports'
str3='Delay between first recover and first case reports'
str2='# Countries reporting'
str4='# Countries reporting'
bx.set_ylabel(str1); dx.set_ylabel(str3)
ax.set_xlabel(r'''Delay to first death and first case from $d_0$''')
cx.set_xlabel(r'''Delay to first recovery and first case from $d_0$''')
bx.set_xlabel(str2); dx.set_xlabel(str2)
ax.set_ylabel(str2); cx.set_ylabel(str4)
ax.legend(loc='upper left')
cx.legend(loc='upper left')
gr.ion(); gr.draw()
fRepDelay.subplots_adjust(left=0.075,bottom=0.075,right=0.9,top=0.95,wspace=0.2,hspace=0.25)
fRepDelay.savefig('./figures_COVID19_dataAnalysis/dam_COVID19_JHU_delays_caseDeaths.png')

#print(aa,bb)


# Now a comparison between those estimates of delays and the delays (not just the first case from the whole world)

# In[17]:


days = np.arange(nDays)
qs = np.arange(0,1,0.01)
delaysFDC = sdFD-sdFC
qDelaysFDC = np.quantile(delaysFDC,qs)
print('Quantiles for delays D-C:',qDelaysFDC)
medianFDCdelay = np.int32(np.median(delaysFDC) )
print('The median delay between first death and first case is %d days'%medianFDCdelay)
print('The mode of the distribution for first death and case delays is %d'% findMode(dFD-dFC))
print('The mode of the distribution for first recovery and case delays is %d'% findMode(dFR-dFC))
print(np.sort(sdFD-sdFC))
print(np.sort(sdFR-sdFC))


# In[18]:


# -----------------------------------
"""Poincaré diagram to plot delays in reports from cumulative counts
of first case to first recovery or death
"""
# -----------------------------------
# prepare some data
v1 = 33; v2=145
delay=18
iStart = np.where(cumC>v1)[0].min()
iStop = np.where(cumC>v2)[0].min()
x=cumC[iStart:iStop]
y=cumD[iStart:iStop]
dd0= days[iStart:iStop]
#
shift=2
gDelay= gr.figure(figsize=(7,5)); gr.ioff()
title="Delay between first case and first death"
gDelay.suptitle(title)
ax = gDelay.add_subplot(111)
ax.plot(days, cumC, 'b', lw=1, alpha=1, label=r'First case report')
ax.plot(days, cumD, 'k', lw=1, alpha=1, label=r'First death report')
d1=delay-shift; d2=delay+shift
ax.plot(days+d1, cumC, 'b', lw=1, alpha=0.5, label=r'First case report + %d days'%d1)
ax.plot(days+d2, cumC, 'b', lw=1, alpha=0.5, label=r'First case report + %d days'%d2)
ax.plot([0,d2+days.max()],[v1,v1],'k:',lw=1)
ax.plot([0,d2+days.max()],[v2,v2],'k:',lw=1)
ax.set_ylabel('# Countries reporting')
ax.set_xlabel('Days from $d_0$')
ax.legend()
gr.ion(); gr.draw();gr.show()
fName="./figures_COVID19_dataAnalysis/dam_COVID19_delaysCumFDC.png"
gDelay.savefig(fName)


# 

# In[19]:


#cFirstCases, binsC = np.histogram(dFC,np.arange(0,(iFD-iFC).max()))
days = np.arange(nDays)
cDelayDeaths2Case, binsD = np.histogram(dFD-dFC,days)
cDelayRecovs2Case, binsR = np.histogram(dFR-dFC,days)
days = np.arange(nDays)

def plotReportingDelays():
    fDelays= gr.figure(figsize=(7,10)); gr.ioff(); rows=2; cols=1;
    ticks = np.arange(0,nDays,7)
    ax=list(); tAx=list()
    fDelays.suptitle('Delays in first case report, first death, and first recovery, all countries')
    for m in range(rows*cols):
        ax.append(fDelays.add_subplot(rows,cols,m+1))
    tAx.append(ax[0].twinx())
    ax[0].bar(binsD[:-1],cDelayDeaths2Case,width=0.8,label='Deaths')
    ax[1].bar(binsR[:-1],cDelayRecovs2Case,width=0.8,label='Recoveries')
    #
    ax[0].set_xlabel('Days first case reported (same country)')
    ax[1].set_xlabel('Days first case reported (same country)')

    for m in range(rows*cols):
        #ax[m].set_yticks(np.arange(0,np.maximum(cDelayRecovs2Case.max()+1,cumC.max()),2))
        ax[m].set_ylabel('# Countries',rotation=90)
        if m>0:
            ax[m].legend(loc='upper right')
        if m<1:
            tAx[m].legend(loc='upper left')
    fDelays.subplots_adjust(left=0.075,bottom=0.075,right=0.9,top=0.95,wspace=0.1,hspace=0.25)
    gr.ion(); gr.draw(); gr.show()
    fDelays.savefig('figures_COVID19_dataAnalysis/dam_COVID19_JHU_delays_AllCountries.png')
    return fDelays

fDelays= plotReportingDelays()


# In[20]:


n=10; print(r'%d countries reported within the first %d days from $d_0$'%(len(sdFC[sdFC<n+1]),n))
countries[iArrival[sdFC<n+1]]


# In[21]:


n=30; print(r'%d countries reported within the first %d days from $d_0$'%(len(sdFC[sdFC<n+1]),n))
countries[iArrival[sdFC<n+1]]


# 

# ## Cases vs deaths in some subsets of countries taking the delay into account

# In[103]:



#
def plotCasesDeathsTS(casesTS,deathsTS,regions,countries, convFactor=1000):
    ii = getIndsRegions(countries, regions)
    figu= gr.figure(figsize=(7,9))
    if convFactor <= 1000:
        figu.suptitle('Deaths vs cases per %d habitants, %s'%(convFactor,dates[-1]))
    elif convFactor == 10**6:
        figu.suptitle('Deaths vs cases per million, %s'%(dates[-1]))
    gr.ioff(); ax=list(); sax=list(); cols=1; rows =3
    ticks= np.arange(0,nDays,7)
    for n in range(len(regions)):
        ax.append(figu.add_subplot(rows,cols,n+1))
        sax.append(inset_axes(parent_axes=ax[n],
                                width="30%", # width = 30% of parent_bbox
                                height="30%", # height : 1 inch
                                loc='lower right'))
        region=ii[n]
        if convFactor <= 1000:
            strDeaths = 'Deaths x {: d}'.format(convFactor)
        elif convFactor == 10**6:
            strDeaths = 'Deaths per million'
        for nn in range(len(region)):
            cas=convFactor*np.float64(casesTS[region[nn]]) /Pops_Millions[countries[region[nn]]]
            dea=convFactor*np.float64(deathsTS[region[nn]]) /Pops_Millions[countries[region[nn]]]
            sax[n].plot(cas,dea,'-',label=countries[region[nn]])
            ax[n].plot(cas,dea,'-',label=countries[region[nn]])
            if convFactor <= 1000:
                ax[n].set_xlabel(r'cases x %d'%convFactor);  ax[n].set_ylabel(r'deaths x %d'%convFactor)
            elif convFactor == 10**6: 
                ax[n].set_xlabel(r'cases per million');  ax[n].set_ylabel(r'deaths per million')
        ymm= ax[n].get_ylim()[1]/3
        xmm= ax[n].get_xlim()[1]/3
        sax[n].set_ylim(0.0,ymm);  sax[n].set_xlim(0.0,xmm);sax[n].set_xticklabels([])
        ax[n].legend(ncol=3,loc='upper left',fontsize=8)
    figu.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.95,wspace=0.2,hspace=0.25)
    gr.ion(); gr.draw(); gr.show()
    return figu


# In[104]:


print(totCases.shape)


# In[105]:


# Setup regions
R1=['China','Japan','Korea, South','Indonesia','Singapore','Australia']
R2=['United Kingdom','Spain','Italy','France','Germany']
#America=['Mexico','US','Argentina','Brazil','Colombia','Chile']
LatinAmerica=['Argentina','Brazil','Colombia','Mexico']
Africa=['Niger','Algeria','Egypt','South Africa']
R3=['US','Canada']
MiddleEast =['Iran','Israel']
#MiddleEast =['Iran','Lebanon', 'West Bank and Gaza','Israel']
regions=[R1,R2+MiddleEast+R3,LatinAmerica+Africa]


# In[112]:


# 
delay = 21
ii = getIndsRegions(countries, regions)
figu=plotCasesDeathsTS(casesTS=totCases[:,:-delay], deathsTS=totDeathCases[:,delay:], regions=regions, countries=countries, convFactor=10**6)
figName='figures_COVID19_dataAnalysis/dam_COVID19_JHU_cases-deaths_x%d_CDdelay%ddays_day%d.png'%(convFactor,delay,nDays)
figu.savefig('./'+figName)


# In[113]:


# 
delay = 14
ii = getIndsRegions(countries, regions)
figu=plotCasesDeathsTS(casesTS=totCases[:,:-delay], deathsTS=totDeathCases[:,delay:], regions=regions, countries=countries, convFactor=10**6)
figName='figures_COVID19_dataAnalysis/dam_COVID19_JHU_cases-deaths_x%d_CDdelay%ddays_day%d.png'%(convFactor,delay,nDays)
figu.savefig('./'+figName)


# In[115]:


# 
delay = 7
ii = getIndsRegions(countries, regions)
figu=plotCasesDeathsTS(casesTS=totCases[:,:-delay], deathsTS=totDeathCases[:,delay:], regions=regions, countries=countries, convFactor=10**6)
figName='figures_COVID19_dataAnalysis/dam_COVID19_JHU_cases-deaths_x%d_CDdelay%ddays_day%d.png'%(convFactor,delay,nDays)
figu.savefig('./'+figName)


# In[114]:


# 
delay = 0
ii = getIndsRegions(countries, regions)
figu=plotCasesDeathsTS(casesTS=totCases, deathsTS=totDeathCases, regions=regions, countries=countries, convFactor=10**6)
figName='figures_COVID19_dataAnalysis/dam_COVID19_JHU_cases-deaths_x%d_CDdelay%ddays_day%d.png'%(convFactor,delay,nDays)
figu.savefig('./'+figName)


# In[22]:


# -------------------
# Search regions to illustrate the case-fatality ratios
# -------------------
Pops_Millions = {'China':1439323776, 'Japan':126476461,'Korea, South':51269185, 'Indonesia':273523615,
                 'Iran':83992949,'Israel':8655535,
                 'Singapore':5850342, 'Mexico':128932753,'US':331002651, 'Canada':37742154, 'Argentina':45195774,
                 'Brazil':212559417, 'Colombia':50882891, 'Niger':24206644, 'Algeria':43851044, 'Egypt':102334404, 
                 'South Africa':59308690,'Spain':46754778,'Italy':60461826,'France':67886011,'Germany':83783942, 
                 'Australia':25499884,'United Kingdom':67886011, 'Venezuela':28870195}


# ## Case-fatality ratios 

# The case fatality ratio is an approximation for the probability of death among cases in an epidemic. In fact, it is an upper bound for the proportion of deaths due to infection, assuming that people that have not been confirmed do not have a higher probability of dying because of the infection.
# 
# The case-fatality ratios can be calculated by dividing each entry in the deaths data frame, by the corresponding entry in the cases data frame.
# 

# In[26]:


def delayCorrectTSPair(a,b,delay=18):
    nPts = len(a)
    a = a[:-DCdelay]
    b = b[DCdelay:]
    return a,b

DCdelay=18; DRdelay=14
DC_totCases,DC_totDeaths=delayCorrectTSPair(totCases.transpose(),totDeathCases.transpose(),delay=DCdelay)
RC_totCases,RC_totRecovs=delayCorrectTSPair(totCases.transpose(),totRecovCases.transpose(),delay=DRdelay)
#
DC_totCases = DC_totCases.transpose()
DC_totDeaths= DC_totDeaths.transpose()
RC_totCases = RC_totCases.transpose()
RC_totRecovs= RC_totRecovs.transpose()


# ### Delayed deaths with respect to case reports in single countries

# In[66]:


def cfr_WholeCountry(totCases,totDeathCases,cou = 'US',insetLims=[60,90],delays=np.arange(5,25,2)):
    iC = np.where(countries == cou)[0][0]
    cas= totCases[iC,:]; dea= totDeathCases[iC,:]
    print(len(cas),len(dea))
    dz = np.zeros(len(dea))
    dz[1:] = dea[1:]-dea[:-1]
    y= np.zeros(len(cas))
    y[1:] = cas[1:]-cas[:-1]

    figCFR1= gr.figure(figsize=(11,9)); gr.ioff(); 
    figCFR1.suptitle('CFRs for different case-death report delays in %s'%cou)
    rows=2; cols=1
    ax=list(); tax=list()
    for n in range(rows*cols):
        ax.append(figCFR1.add_subplot(rows,cols,n+1))
        tax.append(inset_axes(parent_axes=ax[n], width="30%", height="30%", loc='upper right'))

    for d in delays:
        ax[0].plot(days[d:],np.float64(dea[d:])/cas[:-d],ms=1,alpha=0.5,label=r'd=%d'%d)
        tax[0].plot(days[d:],np.float64(dea[d:])/cas[:-d],ms=1,alpha=0.5)    
    for d in delays:
        ax[1].plot(days[d:],np.float64(dz[d:])/y[:-d],ms=1,alpha=0.5,label=r'd=%d'%d)
        tax[1].plot(days[d:],np.float64(dz[d:])/y[:-d],ms=1,alpha=0.5)
    ax[1].set_xlabel(r'Days from $d_0$')
    ax[0].set_ylabel('CFR from cumulative data')
    ax[1].set_ylabel('CFR from incidence')
    for n in range(rows*cols):
        ax[n].legend(loc='upper left',ncol=1)
        tax[n].set_xlim(insetLims)
        ax[n].set_ylim(0,1);tax[n].set_ylim(0.0,0.3)
    figCFR1.subplots_adjust(left=0.075,bottom=0.1,right=0.92,top=0.93,wspace=0.1,hspace=0.1)
    gr.ion(); gr.draw(); gr.show()
    figName='./figures_COVID19_dataAnalysis/dam_COVID19_CFR_delay_%s.png'%(cou)
    figCFR1.savefig(figName)
    return figCFR1


# In[67]:


fCFRCou= cfr_WholeCountry(totCases,totDeathCases,cou = 'Italy',insetLims=[60,nDays],delays=np.arange(5,25,3))


# In[68]:


fCFRCou= cfr_WholeCountry(totCases,totDeathCases,cou = 'Spain',insetLims=[60,nDays],delays=np.arange(5,25,3))


# In[69]:


fCFRCou= cfr_WholeCountry(totCases,totDeathCases,cou = 'US',insetLims=[60,nDays],delays=np.arange(5,25,3))


# In[70]:


fCFRCou= cfr_WholeCountry(totCases,totDeathCases,cou = 'China',insetLims=[40,55],delays=np.arange(5,25,3))


# In[32]:


def zeroCorrectCFR(c,d):
    c[c==0]=1
    print(c.shape,d.shape)
    return d/c

print(totCases.shape,totDeathCases.shape)
print(DC_totCases.shape,DC_totDeaths.shape)
cfr=zeroCorrectCFR(c=np.float64(totCases),d=totDeathCases)
DC_cfr=zeroCorrectCFR(c=np.float64(DC_totCases),d=DC_totDeaths)
print(cfr.shape,DC_cfr.shape)

DC_cfr.max()


# Function to plot the CFR in different subsets of countries chosen specifically to illustrate different dynamics

# In[76]:


#
def plotCFRTS(cfr,dates,regions,countries,move2start=1):
    print('Taking into account %d days'%nDays)
    cfr=100*cfr
    delayDays = np.arange(DCdelay,nDays)
    ii = getIndsRegions(countries, regions)
    figu = gr.figure(figsize=(7,9))
    figu.suptitle('Percentage of dead/confirmed between %s-%s'''%(dates[0],dates[-1]))
    ax=list(); gr.ioff()
    cols=1; rows = 3
    ticks= np.arange(0,nDays,7)
    for n in range(len(regions)):
        ax.append(figu.add_subplot(rows,cols,n+1))
        region=ii[n];
        for nn in range(len(region)):
            thisCFR=cfr[region[nn]]
            ax[n].set_xticks(ticks)
            if move2start:
                startInd = np.maximum( np.where(thisCFR>0)[0].min(),0)
                ax[n].plot(thisCFR[startInd:],'-',label=countries[region[nn]])
                ax[n].set_xlabel('Days from first reported case')
            else:
                ax[n].plot(cfr[region[nn]],'-',label=countries[region[nn]])
                ax[n].set_xticklabels(dates[ticks],{'fontsize':8})
                for label in ax[n].get_xticklabels():
                    label.set_rotation(45)
                    label.set_horizontalalignment('center')
                    label.set_fontsize(8)
        #ax[n].set_xlim(ximin,len(dates))
        ax[n].set_ylim(0,15);
        ax[n].legend(ncol=4,loc='upper left',fontsize=8)
    figu.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.95,wspace=0.2,hspace=0.25)
    gr.ion(); gr.draw(); gr.show()
    strCFR='figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfr_fromFirstLocalCase.png'
    figu.savefig('./'+strCFR)
    return figu


# Search regions to illustrate the case-fatality ratios

# Now let us plot the case fatality ratios of a few countries with reported cases. 

# In[73]:


# ---------------------------------------------
# All relative to the starting days of the pandemia
figu=plotCFRTS(DC_cfr,dates,regions,countries, move2start=1)


# In[35]:


# ---------------------------------------------
# From the day the first case was reported
figu=plotCFRTS(cfr,dates,regions,countries, move2start=0)
strCFR='figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfr_relative2d0.png'
figu.savefig('./'+strCFR)
# ---------------------------------------------


# It is important to consider that the first few reports of deaths usually are biased by the fact that those cases are almost 

# ### Case fatality ratios in detail for countries reporting cases by province
# 
# The countries where there are reports by province are China, United Kingdom, Australia, and Canada. To see this, print the list of countries including repetitions.

# In[36]:


countries_Cases


# The rows that contain the data from China, for instance, are between 42 and 82, inclusive.

# ### CFR analysis for China, UK, and Australia, taking into account the data by province
# 
# The data will be separated into dictionaries, one for each country. 

# China

# In[37]:


China= {'Country name':'China'}
China['cfrWC']=cfr[np.where(countries=='China')[0][0]]
China['cases'], China['indsCases']= gatherDataSingleCountry(cases,'China')
China['deathCases'], China['indsDeaths']= gatherDataSingleCountry(deathCases,'China')
China['cfrs'] = correctedArrayRatio(China['deathCases'],China['cases'])
China['provinces']=cases.iloc[China['indsCases'],0].to_numpy()
China['nProvinces'] = len(China['provinces'])
China['startDaysCases']= findCaseStarts(places=China['provinces'],cases=China['cases'])


# UK

# In[38]:


UK =  {'Country name':'UK'}
UK['cfrWC']=cfr[np.where(countries=='United Kingdom')[0][0]]
UK['cases'], UK['indsCases']= gatherDataSingleCountry(cases,'United Kingdom')
UK['deathCases'], UK['indsDeaths']= gatherDataSingleCountry(deathCases,'United Kingdom')
UK['cfrs'] = correctedArrayRatio(UK['deathCases'],UK['cases'])
UK['provinces']=cases.iloc[UK['indsCases'],0].to_numpy()
UK['nProvinces'] = len(UK['provinces'])
UK['startDaysCases']= findCaseStarts(places=UK['provinces'],cases=UK['cases'])
print(UK['provinces'])
UK['provinces'][6]= 'Great Britain'
print(UK['provinces'])


# Australia

# In[39]:


Australia =  {'Country name':'Australia'}
Australia['cfrWC']=cfr[np.where(countries=='Australia')[0][0]]
Australia['cases'], Australia['indsCases']= gatherDataSingleCountry(cases,'Australia')
Australia['deathCases'], Australia['indsDeaths']= gatherDataSingleCountry(deathCases,'Australia')
Australia['cfrs'] = correctedArrayRatio(Australia['deathCases'],Australia['cases'])
Australia['provinces']=cases.iloc[Australia['indsCases'],0].to_numpy()
Australia['nProvinces'] = len(Australia['provinces'])
Australia['startDaysCases']= findCaseStarts(places=Australia['provinces'],cases=Australia['cases'])


# In[40]:


Australia['startDaysCases']


# Plot comparisons for the three places with averages by region and comparison with the whole country average

# In[41]:


def plotCFRTS_Provinces(place,dates,move2start=1):
    figu = gr.figure(figsize=(7,5))
    figu.suptitle('Percentage of dead/confirmed between %s-%s in %s'''%(dates[0],dates[-1],place['Country name']))
    gr.ioff();cols=1; rows = 1
    ax=figu.add_subplot(rows,cols,1)
    print(place['Country name'])
    ticks= np.arange(0,nDays,7)
    si = place['startDaysCases'].min()
    for nn in range(len(place['provinces'])):
        #print(place['provinces'][nn])
        thisCFR=100*place['cfrs'][nn]
        ax.set_xticks(ticks)
        if move2start==1:
            ax.plot(thisCFR[si:],'-',label=place['provinces'][nn])
            ax.set_xlabel('Days from first reported case')
        else:
            ticks= np.arange(si,nDays,7)
            ax.plot(thisCFR,'-',label=place['provinces'][nn])
            ax.set_xticklabels(dates[ticks],{'fontsize':8})
            for label in ax.get_xticklabels():
                label.set_rotation(0)
                label.set_horizontalalignment('center')
                label.set_fontsize(8)
    avgCFR=100*place['cfrs'].mean(0)
    if move2start==1:
        ax.plot(avgCFR[si:],'k-',alpha=1, lw=3,label='Average CFR from provinces in %s'%place['Country name'])
        ax.plot(100*place['cfrWC'][si:],'k-',alpha=0.35, lw=5,label='CFR from total cases')
        strCFR='figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfr_Provinces'+place['Country name']+'_fromDate0.png'
    else:
        ax.plot(avgCFR,'k-',alpha=1, lw=3,label='Average CFR from provinces in %s'%place['Country name'])
        ax.plot(100*place['cfrWC'],'k-',alpha=0.35, lw=5,label='CFR from total cases')
        strCFR='figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfr_Provinces'+place['Country name']+'_fromFirstLocalReport.png'
    ax.set_ylim(0,10);
    ax.legend(ncol=4,loc='upper center',fontsize=8)
    figu.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.92,wspace=0.2,hspace=0.25)
    figu.savefig('./'+strCFR)
    gr.ion(); gr.draw(); gr.show()
    return figu


# In[42]:


figuChina=plotCFRTS_Provinces(China,dates,move2start=0)


# In[43]:


figuUK=plotCFRTS_Provinces(UK,dates,move2start=0)
print(UK['provinces'])


# In[44]:


figuAust=plotCFRTS_Provinces(Australia,dates,move2start=0)


# ### Predictions of CFR by age groups based on estimates from China and Italy

# The rates of mortality among cases by age group have been reported in some of the countries where COVID-19 has been detected.  The probability that a person with positive SARS-CoV-2 testing dies in the hospital depends on many different factors. Some of these factors depend on the individual, including age, comorbidities like hypertension, obesity and diabetes, etcetera [(10.1016/S0140-6736(20)30566-3)](https://www.sciencedirect.com/science/article/pii/S0140673620305663). Other factors hinge on the availability of health care, and the quality of care in different countries.  Of note, the mortality in Italy, South Korea, and China, are constrasting cases of study of the known patterns patterns for CFRs results due to different conditions in these countries, different population structures, and health care policies. 
# 
# As a means to generate estimates of the CFRs from different age groups, we propose a simple calculation of proportions based on cases. Depending on an unknown factor that would account for under reporting, these estimates can be used to generate time series of CFRs by age. 

# The age groups to take into account are \[0,10),\[10,20), \[20,30), \[30,40), \[40,50), \[50,60), \[60,70), \[70,80), 80+. 
# 
# The data that will be used was last estimated in February 11, 2020 for China, March 17, 2020, for Italy, and March 24, 2020 for South Korea [(video of Kim Woo-Joo from the Guro University Hospital Guro, University of Korea (March 27, 2020)](https://www.youtube.com/watch?v=gAk7aX5hksU&feature=youtu.be).
# 
# There are reports of up to 20% of symptomatic cases, and in South Korea, for instance, about 50% of hospitalizations are people 60 years old or older, it is important to obtain numbers to estimate the size of a population in need for critical care. 

# In[45]:


SKorea={'ageCFR': np.array([0,0,0,0.001,0.001,0.004,0.015,0.063,0.116]),'CFR_20200324':0.02,'Country name':'Korea, South'}
Italy={'ageCFR': np.array([0,0,0,0.003,0.004,0.01,0.035,0.128,0.202]),'CFR_20200317':0.072,'Country name':'Italy'}
China['ageCFR'] = np.array([0,0.002,0.002,0.002,0.004,0.013,0.036,0.08,0.14])
China['CFR_20200211'] = 0.02
aGroups=np.arange(0,90,10)


# In[46]:


def sigmoid(a, aMax=0.1,a0=60.0,n=2): 
    aa = a**n
    return aMax* aa /(aa + a0**n)


# In[47]:


ageCountries=[China, SKorea, Italy]
ages = np.arange(0,90)
f11=gr.figure(figsize=(9,5)); gr.ioff(); rows =1;cols=1
f11.suptitle('CFR by age')
lshift=2.5;  ax=f11.add_subplot(rows,cols,1)
for n in range(len(ageCountries)):
    ax.bar(aGroups+ (n*lshift), 100*ageCountries[n]['ageCFR'],width=2,align='edge',label=ageCountries[n]['Country name']);
ax.plot(ages, 100*sigmoid(a=ages,a0=75,aMax=0.22,n=9))
ax.plot(ages, 100*sigmoid(a=ages,a0=80,aMax=0.17,n=9))
ax.plot(ages, 100*sigmoid(a=ages,a0=82,aMax=0.33,n=9))
ax.legend(loc='upper left')
ax.set_xticks(aGroups)
ax.set_yticks(np.arange(0,20+1,2))
ax.set_xlabel('Age groups')
gr.ion(); gr.draw();  gr.show()
f11Name='figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfrByAge_China+SKorea+Italy.png'
f11.savefig(f11Name)


# These profiles are similar in that they are all increasing as a function of age, reaching a half-maximum height at ages between 50 and 60 years of age. 

# We can establish qualitative estimations of the percentage of deaths based on these age CFRs. To do so, we calculate the weights of the CFRs in each age group, which represent the proportion of deaths in each group.
# 
# #### Calculation of proportions of death within cases, by age groups
# 
# Let $c_i$ represent the case-fatality ratio in $i$th age group, $i \in \left\{0,1,...,8\right\}$
# The weight of the $i$th age group on the total death reported via CFRs is
# $$w_i = \frac{c_i}{\sum_{j=0}^{n-1}} w_j.$$
# Then the probability of death in the $i$th age group can be estimated by multiplying the CFR and the weight of the $i$th age group on the CFR. Then, multiplication by the number of cases would yield the death toll by age group. 
#  

# A few sutile differences can be noted from these data. For instance, the CFRs for the first 3 groups in Italy and South Korea are very similar. In Italy, the percentage of dead cases is very large among the last two age groups, which in part is due to the fact that Italians implemented a triage system that sometimes denied respirators to patients with small chances of survival. 

# In[48]:


ageCountries=[China,SKorea,Italy]
for n in range(len(ageCountries)):
    C = ageCountries[n]
    print(C['Country name'])
    ageCFRweights = C['ageCFR'].sum()
    C['ageDeadCaseProps'] = C['ageCFR']/ageCFRweights
    C['cfrTotCases']=cfr[np.where(countries==C['Country name'])[0][0]]
    C['totCases'] = totCases[np.where(countries==C['Country name'])[0][0]]
    C['totDeathCases'] = totDeathCases[np.where(countries==C['Country name'])[0][0]]
    C['totDeathCases_ageProps'] = np.zeros((len(aGroups),len(C['totDeathCases'])),'int64')
    for nn in range(len(aGroups)):
        C['totDeathCases_ageProps'][nn,:]= C['ageDeadCaseProps'][nn] * C['totCases']* C['cfrTotCases']
        print(C['totDeathCases_ageProps'][nn])



# In[49]:


ages = np.arange(0,90)
f10b=gr.figure(figsize=(7,9)); gr.ioff(); rows =2;cols=1
f10b.suptitle('Age related CFRs and their relative weights'); ax=list()
for m in range(rows*cols):
    ax.append(f10b.add_subplot(rows,cols,m+1))
lshift=2.5;  
for n in range(len(ageCountries)):
    C=ageCountries[n]
    ax[0].bar(aGroups+ (n*lshift), 100*C['ageCFR'],width=2,align='edge',label=ageCountries[n]['Country name']);
    ax[1].bar(aGroups+ (n*lshift), 100*C['ageDeadCaseProps'],width=2,align='edge',label=ageCountries[n]['Country name']);
ax[0].plot(ages, 100*sigmoid(a=ages,a0=75,aMax=0.22,n=9))
ax[0].plot(ages, 100*sigmoid(a=ages,a0=79,aMax=0.17,n=9))
ax[0].plot(ages, 100*sigmoid(a=ages,a0=81,aMax=0.33,n=9))
ax[1].plot(ages, 100*sigmoid(a=ages,a0=75,aMax=0.78,n=9))
ax[1].plot(ages, 100*sigmoid(a=ages,a0=80,aMax=0.9,n=9))
ax[1].plot(ages, 100*sigmoid(a=ages,a0=84,aMax=0.92,n=9))
for m in range(rows*cols):
    ax[m].legend(loc='lower left'); 
    ax[m].set_xticks(aGroups)
    ax[m].set_xlabel('Age groups')
ax[0].set_yticks(np.arange(0,26,5))
ax[1].set_yticks(np.arange(0,66,5))
ax[0].text(0,20,'CFRs by age (percent)')
ax[1].text(0,60,'Contributions to death within cases by age (percent)')
f10b.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.9,wspace=0.2,hspace=0.25)
gr.ion(); gr.draw();  gr.show()
f10bName='figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfr+propDeathCases_ByAge_China+SKorea+Italy_OneFigure.png'
f10b.savefig(f10bName)


# In[50]:


ages = np.arange(0,90)
f12=gr.figure(figsize=(13,9)); gr.ioff(); rows =3;cols=1
f12.suptitle('Deaths due to COVID-19 by age')
ax=list()
ticks= np.arange(0,nDays,3)
convFactor =1
for n in range(len(ageCountries)):
    ax.append(f12.add_subplot(rows,cols,n+1))
    C = ageCountries[n]
    print(C['Country name'])
    ax[n].plot(C['totDeathCases'],'k--', alpha=1, lw=1,label='Deaths reported '+ C['Country name']);
    ax[n].plot(C['totDeathCases_ageProps'].sum(0),'k',lw=1,label='Sum of estimates '+ C['Country name']);
    for nn in range(len(aGroups)):
        ax[n].plot(C['totDeathCases_ageProps'][nn],label=C['Country name']+ '[%d,%d)]'%(aGroups[nn],aGroups[nn]+10));
    ax[n].legend(loc='upper left',ncol=2, fontsize=8)
    ax[n].set_xticks(ticks)
    ax[n].set_xticklabels(dates[ticks],{'fontsize':8})
    for label in ax[n].get_xticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment('center')
        label.set_fontsize(8)
f12.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.95,wspace=0.2,hspace=0.25)
gr.ion(); gr.draw();  gr.show()

f12Name='figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfr+propDeathCasesByAgeTS.png'
f12.savefig(f12Name)

print('%d days'%nDays)


# In[51]:


print(SKorea['totDeathCases_ageProps'].sum(0)[-1])
print(Italy['totDeathCases_ageProps'].sum(0)[-1])
print(China['totDeathCases_ageProps'].sum(0)[-1])


# ### Projections for Mexico

# In[52]:


Mexico = dict()
Mexico['cfr'] =cfr[np.where(countries=='Mexico')[0][0]]
Mexico['totCases']=totCases[np.where(countries=='Mexico')[0][0]]
Mexico['totDeathCases']=totDeathCases[np.where(countries=='Mexico')[0][0]]
Mexico['startInd'] = np.where(Mexico['totCases']>0)[0].min()
Mexico['popMillions']= Pops_Millions['Mexico']
Mexico['delayFromReport0']=np.where(Mexico['totCases']>0)[0].min()
print('First case reported %d days after Report 0'%Mexico['delayFromReport0'])


# In[53]:


def estimateDeathsByAgeMexico(Mexico, dates, subReportFactor=1):
    si = Mexico['delayFromReport0']
    ticks= np.arange(0,len(dates),7)
    casesMexico = subReportFactor * Mexico['totCases']
    f13=gr.figure(figsize=(7,9)); gr.ioff(); rows =3;cols=1;ax=list()
    for n in range(len(ageCountries)):
        ax.append(f13.add_subplot(rows,cols,n+1))
        C = ageCountries[n]
        ax[n].plot(subReportFactor* Mexico['totDeathCases'],'k',lw=1,label='Deaths reported in Mexico');
        ax[n].set_title('Projection based on contributions by age groups from %s using a subreport factor of %d'%(C['Country name'],subReportFactor))
        for nn in range(len(aGroups)):
            ff = casesMexico * Mexico['cfr'] * C['ageDeadCaseProps'][nn]
            ax[n].plot(ff,label='Mexico [%d,%d)]'%(aGroups[nn],aGroups[nn]+10));
            ax[n].legend(loc='upper left',ncol=3, fontsize=8)
        ax[n].set_xticks(ticks)
        ax[n].set_xticklabels(dates[ticks],{'fontsize':8})
        for label in ax[n].get_xticklabels():
            label.set_rotation(0)
            label.set_horizontalalignment('center')
            label.set_fontsize(8)
        ax[n].set_xlim(si,len(Mexico['totDeathCases']))
    f13.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.95,wspace=0.1,hspace=0.25)
    f13Name='figures_COVID19_dataAnalysis/dam_COVID19_JHU_cfr+propDeathCasesByAgeTS_EstimatesMexico_subReportFactor%d.png'%subReportFactor
    gr.ion(); gr.draw();  gr.show()
    f13.savefig(f13Name)
    return f13



# ### Estimates with a conversion factor of 1

# In[54]:


f13_1 = estimateDeathsByAgeMexico(Mexico, dates,subReportFactor=1)


# ### Estimates with a conversion factor of 10

# In[55]:


f13_10 = estimateDeathsByAgeMexico(Mexico, dates,subReportFactor=10)


# ### Estimates with a conversion factor of 12

# In[56]:


f13_12 = estimateDeathsByAgeMexico(Mexico, dates,subReportFactor=12)


# In[ ]:





# ## Analysis of the dynamics assuming SIR-like evolution
# 
# Assume that the population is divided in three subsets representing the non-infected and susceptible, the infected, and those that can no-longer participate in the chain of infections due to immunity or death. Let the densities of those subpopulations be represented by $x$, $y$, and $z$ respectively, with $1=x+y+z$. The dynamics for $(x,y,z)$ can be assume to follow 
# \begin{eqnarray}
# \partial_t y &=& \beta y \bigg( 1- z - y \bigg) - y \bigg(\gamma + \delta \bigg)\\
# \partial_t z &=& \left( \gamma + \delta \right)y
# \end{eqnarray}
# 
# 
# If $h$ is the time step for sampling, let $y_{n}$ represent the number of cases at time $nh$. From data, the cumulative number of cases at step $nh$ can be written as 
# $$C_{n} = \sum_{k=0}^{n} y_k$$
# which means that the number of new cases (incidence) at time nh is 
# $$y_n = C_n - C_{n-1}.$$
# Similarly, for the cumulative recoveries and deaths, respectively
# $$R_{n}+D_{n} = \sum_{k=0}^{n} z_k,$$
# so 
# $$z_n = (R_{n}-R_{n-1})+(D_{n}-D_{n-1}).$$
# 
# From the SIR classical formulation, assuming $\gamma$ is the rate of removal by recovery and acquisition of immunity and $\delta$ is the fatality rate due to infection, the time dependent change in the prevalence is given by
# \begin{eqnarray}
# \partial_t y &=& y  \beta \big( 1- z - y  \big) - \partial_t z
# \end{eqnarray}
# with $z$ representing the sum of the recoveries and deaths at time $t= kh$.
# The change in $y$ can then be written as 
# $$\partial_t y =  \beta y \bigg( \alpha -z - y  \bigg), $$
# with $$\alpha = 1 -\frac{\gamma + \delta}{\beta}.$$
# This means that the incidence and the prevalence are related cuadratically in the $(y,\partial_t y)$ plane.
# Using the data, the change in $y$ can be written as 
# $$\frac{y_{k+1}-y_{k}}{h} \approx \beta y \bigg( \alpha -z - y \bigg),$$
# which means that 
# $$y_{k+1}-y_{k} \approx h \beta y \bigg( \alpha -z - y \bigg).$$
# 
# #### Parameters for modeling 
# From the equation for $z$,
# $$\gamma + \delta = \frac{\partial_t z}{y} \approx \frac{z_{k+1}-z_k}{h y_k}$$ 
# Also, 
# $$\beta = \frac{\partial_t y + \partial_t z}{y \bigg( 1 -z - y  \bigg)} $$ 
# 
# 

# #### Parabola describing the incidence as a function of the cases

# In[59]:


#def casePhasePlane(b=50.0, cou='Italy'):
def casePhasePlane(cou='Italy',a=0.005,b=10):
    cInd = np.where(countries==cou)[0][0]
    tCases= np.float64(totCases[cInd,:])/Pops_Millions[cou]
    tDeaths=np.float64(totDeathCases[cInd,:])/Pops_Millions[cou]
    tRecovs=np.float64(totRecovCases[cInd,:])/Pops_Millions[cou]
    y = tCases - tRecovs - tDeaths
    z = tRecovs + tDeaths
    dz=np.zeros(len(z))
    dz[1:]= z[1:]-z[:1]
    dy = np.zeros(len(y))
    dy[1:]= y[1:]-y[:-1]
    gammaDelta = dz/y
    beta = (dy+dz)/((1-z-y)*y)
    a = a+z
    fy = b*y*(a-z-y) 
    #
    ff =gr.figure(figsize=(17,5)); gr.ioff(); rows=1; cols=3
    ax1= ff.add_subplot(rows,cols,1); tax1= ax1.twinx()
    ax2= ff.add_subplot(rows,cols,2); 
    ax3= ff.add_subplot(rows,cols,3); 
    ax1.plot(days, y,'o',markeredgecolor='blue', markerfacecolor='white',label=r'$y_k$')
    #ax1.plot(days, tCases,'.',label=r'$(k,C_k)$')
    tax1.plot(days, beta/gammaDelta,'.',label=r'$\beta/\gamma$')
    ax1.plot(days, tDeaths,'.',label=r'$D_k$')
    ax1.plot(days, tRecovs,'.',label=r'$R_k$')
    ax1.plot(days, z,'.',label=r'$z_k=D_k + R_k$')
    ax2.plot(dy,y,'o',markeredgecolor='blue', markerfacecolor='white',label=r'$(y_k - y_{k-1}, y)$')
    ax2.plot(fy,y,'.',color='orange',label=r'$(\beta (a-z-y)y,y)$')
    ax3.plot(1-z-y,y,'.',label=r'$(1-z-y, y)$')
    ax2.plot([0,0],[0,y.max()],'k:')
    ax1.set_ylim(0,1.1* y.max())
    ax2.set_ylim(0,1.1* y.max())    
    ax2.set_xlim(dy.min(),1.3*dy.max())
    ax1.set_xlabel(r'$t$ (days)');ax2.set_xlabel(r'$y_k - y_{k-1}$')
    ax1.set_ylabel(r'$y_k$');ax2.set_ylabel(r'$y_k$')
    ax2.legend(loc='lower right'); ax3.legend(loc='upper right'); 
    ax1.legend(loc='center left'); tax1.legend(loc='upper left')
    gr.ion(); gr.draw()
    ff.suptitle('Phase plane, %s'%cou)
    ff.subplots_adjust(left=0.075,bottom=0.075,right=0.9,top=0.92,wspace=0.2,hspace=0.25)
    ff.savefig('figures_COVID19_dataAnalysis/dam_COVID19_JHU_phasePlane+%s.png'%cou)
    return y, z, dy, dz, beta, gammaDelta, ff

def gammaDelta_Estimate(y, dz):
    return (dz/y)

def betaEstimate(y,z,dy, dz):
    return (dy+dz)/((1-z-y)*y)


# In[60]:


y,z, dy,dz,beta,gammaDelta,ff=casePhasePlane(cou='Italy',a=0.0018,b=90)


# In[61]:


y,z, dy,dz,beta,gammaDelta,ff=casePhasePlane(cou='Spain',a=0.002,b=110)


# In[62]:


y,z, dy,dz,beta,gammaDelta,ff=casePhasePlane(cou='US',a=0.0023,b=70)


# In[63]:


#y,z, dy,dz,beta,gammaDelta,ff=casePhasePlane(cou='Italy',a=0.00175,b=90)
y,z, dy,dz,beta,gammaDelta,ff=casePhasePlane(cou='China',a=0.00004,b=6000)
#print((1-y-z)*beta/gammaDelta )


# In[64]:


y,z, dy,dz,beta,gammaDelta,ff=casePhasePlane(cou='Mexico',a=0.00006,b=4000)


# In[ ]:




