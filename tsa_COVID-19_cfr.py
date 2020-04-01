# COVID-19_cfr_Jan-March_2020

import pandas as pd
import numpy as np
import matplotlib.pylab as gr
from matplotlib.dates import mdates
import datetime

ref="""Data obtained from https://github.com/CSSEGISandData/COVID-19/"""

urlData=1
if urlData==1:
    srcDir='https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/'
    urlCases = srcDir+'time_series_covid19_confirmed_global.csv'
    cases = pd.read_csv(urlCases,index_col=None)
    urlDeaths = srcDir+'time_series_covid19_deaths_global.csv'
    deaths = pd.read_csv(urlDeaths,index_col=None)
    urlRecov = srcDir+'time_series_covid19_recovered_global.csv'
    recov = pd.read_csv(urlRecov,index_col=None)

# From local directory
localData=0
if localData==1:
    srcDir='./COVID-19/csse_covid_19_data/csse_covid_19_time_series/'
    casesF='time_series_covid19_confirmed_global.csv'
    deathsF='time_series_covid19_deaths_global.csv'
    recovF='time_series_covid19_recovered_global.csv'
    cases=pd.read_csv(srcDir+casesF)
    deaths=pd.read_csv(srcDir+deathsF)
    recov=pd.read_csv(srcDir+recovF)
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
#base = datetime.datetime(2020, 1, 22)
#dates = np.array([base + datetime.timedelta(days=i) for i in range(nDays)])
dates = cases.columns[4:]
nDays = len(dates)
print('Dates:\n',dates)

# -------------------
print("""Data curation
- Some countries report by province.
- Arrays that include whole country cases, deaths, and recovered obtained
and separated for analysis
- Recovered are treated separately because reports do not include as many places as cases or deaths
- Case fatality ratios are calculated to analyse the epidemia in different provinces and contrast with whole country estimates
""")
# -------------------
npCases = cases.to_numpy()
countries = np.unique(npCases[:,1])
nCountries = len(countries)
# -------------------
# Ordering data to gather same country data
# -------------------
def sortDataByCountry(df,nHeaderCols=4):
    iCoun= df.iloc[:,1].to_numpy().argsort()
    x = df.iloc[iCoun,nHeaderCols:].to_numpy()
    places = df.iloc[iCoun,1].to_numpy()
    return x, places

gCases, countries_Cases= sortDataByCountry(cases, nHeaderCols)
gDeaths,countries_Deaths = sortDataByCountry(deaths, nHeaderCols)
gRecov,countries_Recov = sortDataByCountry(recov, nHeaderCols)

# -------------------
# Sum the counts from each country and construct a new array
# -------------------
def gatherDataByCountry(df,nHeaderCols=4):
    cc = df.iloc[:,1].to_numpy()
    countries = np.unique(df.iloc[:,1].to_numpy())
    nCountries = len(countries);
    x = list()
    for n in range(nCountries):
        iC= np.where(cc== countries[n])[0]
        a =df.iloc[iC,nHeaderCols:].to_numpy()
        x.append(a.sum(0))
    return np.array(x)

# These arrays have the same size as the countries array (unique countries)
totCases=gatherDataByCountry(df=cases,nHeaderCols=4)
totDeaths=gatherDataByCountry(df=deaths,nHeaderCols=4)
totRecov=gatherDataByCountry(df=recov,nHeaderCols=4)

# -------------------
# Ratios by array. Correction in cases there are zeros in the denominators
# -------------------
def correctedArrayRatio(a,b):
    bCorr= b.copy()
    bCorr[bCorr==0]=1
    return a/np.maximum(b,bCorr)


# -------------------
# Case-Fatality ratio analysis.
# Total cases
# -------------------
cfrAnalysis="""The case fatality ratio is an approximation for the probability of death among cases in an epidemic.
In fact, it is an upper bound for the proportion of deaths due to infection, assuming that people that have not been confirmed do not have a higher probability of dying because of the infection"""
print(cfrAnalysis)
"""
cfr= correctedArrayRatio(totDeaths,totCases)
# Search regions to illustrate the case-fatality ratios
Asia=['China','Japan','Indonesia','Vietnam']
Oceania=['Australia','New Zealand']
Northamerica=['Mexico','US','Canada']
Southamerica=['Argentina','Brazil','Colombia','Guatemala']
Africa=['Niger','Algeria','South Africa']
Europe=['Spain','Italy','Germany','United Kingdom']
regions=[Asia,Oceania,Northamerica,Southamerica,Africa,Europe]
# Indexing and gathering of data

def getIndsRegion(region):
    i = [np.where(countries==region[nn])[0][0] for nn in range(len(region))]
    return i

def getIndsRegions(regions):
    ii =list()
    for mm in range(len(regions)):
        reg = regions[mm]
        ii.append(getIndsRegion(reg))
    return ii

ii = getIndsRegions(regions)

for n in range(rows):
    j = ii[n]
    for nn in range(len(ii[n])):
        print(cfr[j[nn]])

figu = gr.figure(figsize=(15,9))
ax=list(); gr.ioff()
rows=len(regions)
ticks=np.arange(0,nDays,7)
for n in range(rows):
    ax.append(figu.add_subplot(rows,1,n+1))
    j=ii[n]
    for nn in range(len(j)):
        ax[n].plot(cfr[j[nn]],'-',label=countries[j[nn]])
        ax[n].plot(cfr[j[nn]],'wo',ms=2,mec='k',mfc='w')
    ax[n].legend(ncol=5,loc='upper left',fontsize=8)
    ax[n].set_xticks(ticks)
    ax[n].set_xticklabels(dates[ticks],{'fontsize':8})
    for label in ax[n].get_xticklabels():
        label.set_rotation(0)
        label.set_horizontalalignment('center')
        label.set_fontsize(8)
figu.subplots_adjust(left=0.1,bottom=0.1,right=0.97,top=0.97,wspace=0.2,hspace=0.25)
gr.ion(); gr.draw(); gr.show()


# ---------------------------------------------
# ---------------------------------------------
# ---------------------------------------------






uCountries = np.unique(countries)
nCountries= len(uCountries)
tCases= list()
tDeaths= list()
#tRecov= list()
for n in range(nCountries):
    totCases,totDeaths= getCountryTotals(uCountries[n],countries,c,d)
    tCases.append(totCases)
    tDeaths.append(totDeaths)
    #tRecov.append(totRecov)

tCases=np.array(tCases)
tDeaths=np.array(tDeaths)

def correctedArrayRatio(a,b):
    bCorr= b.copy()
    bCorr[bCorr==0]=1
    return a/np.maximum(b,bCorr)

def totalCFRCRR(country,countries):
    iCountry = np.where(countries==country)[0]
    tc = c[iCountry].sum(0)
    td = d[iCountry].sum(0)
    tr = r[iCountry].sum(0)
    cfr = correctedArrayRatio(totDeaths,totCases)
    crr = correctedArrayRatio(totRecov,totCases)
    tot = {'cases':tc,'deaths':td,'decov':tr,
    'cfr':cfr,'crr':crr,'indices':iCountry}
    return tot



# -----------------------------------
# China analysis only
# -----------------------------------
totChina= getCountryTotals(countries,'China')
iChina = np.where(countries=='China')[0]

f=gr.figure(figsize=(17,9)); gr.ioff()
for n in iChina:
    gr.plot( 100*cfr[n],lw=1, alpha=0.5, label=provinces[n])
gr.plot(100*cfr.mean(0),'r', alpha=0.8, lw=2, label='avg provinces')
gr.plot(100*totChina['cfr'],'k', lw=3, alpha=0.8, label='All China')
gr.legend(ncol=4,fontsize=8)
gr.ylim(0,5)
gr.xlabel('Days from January 22')
gr.ylabel('% death/conf ')
gr.ion(); gr.draw(); gr.show()
