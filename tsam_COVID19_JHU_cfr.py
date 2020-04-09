# tsam_COVID19_cfr_JHU
import os
cwd = os.getcwd()
from tsam_COVID19_baseCode import *

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

if 1:
    urlTests='https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/testing/'
    testsF='covid-testing-all-observations.csv'
    tests = pd.read_csv(urlTests+testsF,index_col=None)

#world population Data
if 1:
    urlWorldPop='https://stats.oecd.org/sdmx-json/data/DP_LIVE/.POP.../OECD?contentType=csv&detail=code&separator=comma&csv-lang=en'
    worldPops = pd.read_csv(urlWorldPop,index_col=None)


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
dates = cases.columns[4:]
nDays = len(dates)
print('Got data from %d Dates:\n'%(len(dates)))
print(dates)

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
# Gather same country data
# -------------------
gCases, countries_Cases= sortDataByCountry(cases, nHeaderCols)
gDeaths,countries_Deaths = sortDataByCountry(deaths, nHeaderCols)
gRecov,countries_Recov = sortDataByCountry(recov, nHeaderCols)

# -------------------
# Sum the counts from each country and construct a new array
# -------------------
# These arrays have the same size as the countries array (unique countries)
totCases=gatherDataByCountry(df=cases,nHeaderCols=4)
totDeaths=gatherDataByCountry(df=deaths,nHeaderCols=4)
totRecov=gatherDataByCountry(df=recov,nHeaderCols=4)

# -------------------
# Search regions to illustrate the case-fatality ratios
# -------------------
Pops_Millions = {'China':1439323776, 'Japan':126476461,'Korea, South':51269185, 'Indonesia':273523615, 'Singapore':5850342, 'Mexico':128932753,
'US':331002651, 'Canada':37742154, 'Argentina':45195774, 'Brazil':212559417, 'Colombia':50882891, 'Niger':24206644, 'Algeria':43851044, 'Egypt':102334404, 'South Africa':59308690,'Spain':46754778,'Italy':60461826,'France':67886011,'Germany':83783942, 'Australia':25499884,'United Kingdom':67886011,'Iran':83992949,'Israel':8655535}
#
R1=['China','Japan','Korea, South','Indonesia','Singapore','Australia']
R2=['United Kingdom','Spain','Italy','France']
#America=['Mexico','US','Argentina','Brazil','Colombia','Chile']
LatinAmerica=['Mexico','Argentina','Brazil','Colombia']
Africa=['Niger','Algeria','Egypt','South Africa']
R3=['US','Canada','Germany']
R4 =['Iran','Israel']
#MiddleEast =['Iran','Lebanon', 'West Bank and Gaza','Israel']
regions=[R1,LatinAmerica+Africa,R3+R4,R2]
#
ii = getIndsRegions(countries, regions)
#
def plotCasesDeathsTS(casesTS,deathsTS,regions,countries, convFactor=1000):
    ii = getIndsRegions(countries, regions)
    figu= gr.figure(figsize=(15,9))
    figu.suptitle('''Deaths vs cases per %d habitants between %s-%s'''%(convFactor,dates[0],dates[-1]))
    gr.ioff(); ax=list(); sax=list(); cols=2
    rows = np.int32(np.ceil(len(regions)/cols))
    ticks= np.arange(0,nDays,7)
    for n in range(len(regions)):
        ax.append(figu.add_subplot(rows,cols,n+1))
        sax.append(inset_axes(parent_axes=ax[n],
                                width="30%", # width = 30% of parent_bbox
                                height="30%", # height : 1 inch
                                loc='lower right'))
        region=ii[n]
        for nn in range(len(region)):
            cas=convFactor*np.float64(casesTS[region[nn]])/Pops_Millions[countries[region[nn]]]
            dea=convFactor*np.float64(deathsTS[region[nn]])/Pops_Millions[countries[region[nn]]]
            ax[n].plot(cas,dea,'-',label=countries[region[nn]])
            ax[n].set_xlabel(r'cases per %d habitants'%convFactor)
            ax[n].set_ylabel(r'deaths per %d habitants'%convFactor)
            sax[n].plot(cas,dea,'-',label=countries[region[nn]])
        ymm= ax[n].get_ylim()[1]/3
        xmm= ax[n].get_xlim()[1]/2
        sax[n].set_ylim(0.0,ymm);  sax[n].set_xlim(0.0,xmm); sax[n].set_xticklabels([])
        ax[n].legend(ncol=6,loc='upper left',fontsize=8)
    figu.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.95,wspace=0.2,hspace=0.25)
    gr.ion(); gr.draw(); gr.show()
    return figu

figu=plotCasesDeathsTS(casesTS=totCases,deathsTS=totDeaths,regions=regions,countries=countries)
strCasesDeaths='tsam_COVID19_cases-deaths1000_JHU.png'
figu.savefig('./'+strCasesDeaths)

# -------------------
# Case-Fatality ratio analysis.
# Total cases
# -------------------
cfrAnalysis= '''The case fatality ratio is an approximation for the probability of death among cases in an epidemic.
In fact, it is an upper bound for the proportion of deaths due to infection, assuming that people that have not been confirmed do not have a higher probability of dying because of the infection'''
print(cfrAnalysis)
#
def plotCFRTS(cfr,dates,regions,countries, move2start=1):
    cfr=100*cfr
    ii = getIndsRegions(countries, regions)
    figu = gr.figure(figsize=(15,9))
    figu.suptitle('Percentage of dead/confirmed between %s-%s'''%(dates[0],dates[-1]))
    ax=list(); gr.ioff()
    cols=2
    rows = np.int32(np.ceil(len(regions)/cols))
    ticks= np.arange(0,nDays,7)
    for n in range(len(regions)):
        ax.append(figu.add_subplot(rows,cols,n+1))
        region=ii[n];
        ax[n].set_xticks(ticks)
        for nn in range(len(region)):
            thisCFR=cfr[region[nn]]
            if move2start:
                startInd = np.maximum( np.where(thisCFR>0)[0].min(),0)
                ax[n].plot(thisCFR[startInd:],'-',label=countries[region[nn]])
                ax[n].set_xlabel('Days from first reported case')
            else:
                ax[n].plot(cfr[region[nn]],'-',label=countries[region[nn]])
                ax[n].set_xticklabels(dates[ticks],{'fontsize':8})
                for label in ax[n].get_xticklabels():
                    label.set_rotation(0)
                    label.set_horizontalalignment('center')
                    label.set_fontsize(8)
        #ax[n].set_xlim(ximin,len(dates))
        ax[n].set_ylim(0,15);
        ax[n].legend(ncol=5,loc='upper left',fontsize=8)
    figu.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.95,wspace=0.2,hspace=0.25)
    gr.ion(); gr.draw(); gr.show()
    return figu
#
cfr= correctedArrayRatio(totDeaths,totCases)
regions=[R1,LatinAmerica+Africa,R3+R4,R2]
# ---------------------------------------------
# All relative to the starting days of the pandemia
figu=plotCFRTS(cfr,dates,regions,countries, move2start=0)
strCFR='tsam_COVID19_cfr_JHU_fromFirstCaseInChina.png'
figu.savefig('./'+strCFR)
# ---------------------------------------------
# From the day the first case was reported
figu=plotCFRTS(cfr,dates,regions,countries, move2start=1)
strCFR='tsam_COVID19_cfr_JHU_fromFirstLocalReport.png'
figu.savefig('./'+strCFR)
# ---------------------------------------------


# ---------------------------------------------
# Analysis of CFRs within countries where there are reports by provinces
# UK, China, Australia
# ---------------------------------------------

# China
China= {'Country name':'China'}
China['cfrTotCases']=cfr[np.where(countries=='China')[0][0]]
China['cases'], China['indsCases']= gatherDataSingleCountry(cases,'China')
China['deaths'], China['indsDeaths']= gatherDataSingleCountry(deaths,'China')
China['cfr'] = correctedArrayRatio(China['deaths'],China['cases'])
China['provinces']=cases.iloc[China['indsCases'],0].to_numpy()
China['nProvinces'] = len(China['provinces'])
China['startDaysCases']= findCaseStarts(places=China['provinces'],cases=China['cases'])
#
UK =  {'Country name':'UK'}
UK['cfrTotCases']=cfr[np.where(countries=='United Kingdom')[0][0]]
UK['cases'], UK['indsCases']= gatherDataSingleCountry(cases,'United Kingdom')
UK['deaths'], UK['indsDeaths']= gatherDataSingleCountry(deaths,'United Kingdom')
UK['cfr'] = correctedArrayRatio(UK['deaths'],UK['cases'])
UK['provinces']=cases.iloc[UK['indsCases'],0].to_numpy()
UK['nProvinces'] = len(UK['provinces'])
UK['startDaysCases']= findCaseStarts(places=UK['provinces'],cases=UK['cases'])
#
Australia =  {'Country name':'Australia'}
Australia['cfrTotCases']=cfr[np.where(countries=='Australia')[0][0]]
Australia['cases'], Australia['indsCases']= gatherDataSingleCountry(cases,'Australia')
Australia['deaths'], Australia['indsDeaths']= gatherDataSingleCountry(deaths,'Australia')
Australia['cfr'] = correctedArrayRatio(Australia['deaths'],Australia['cases'])
Australia['provinces']=cases.iloc[Australia['indsCases'],0].to_numpy()
Australia['nProvinces'] = len(Australia['provinces'])
Australia['startDaysCases']= findCaseStarts(places=Australia['provinces'],cases=Australia['cases'])
#
def plotCFRTS_Provinces(places,dates,move2start=1):
    figu = gr.figure(figsize=(15,9))
    figu.suptitle('Percentage of dead/confirmed between %s-%s'''%(dates[0],dates[-1]))
    ax=list(); gr.ioff()
    cols=1
    rows = len(places)
    for n in range(rows*cols):
        ax.append(figu.add_subplot(rows,cols,n+1))
        place = places[n]
        ticks= np.arange(0,nDays,7)
        si = place['startDaysCases'].min()
        for nn in range(len(place['provinces'])):
            thisCFR=100*place['cfr'][nn]
            ax[n].set_xticks(ticks)
            if move2start==1:
                ax[n].plot(thisCFR[si:],'-',label=place['provinces'][nn])
                ax[n].set_xlabel('Days from first reported case')
            else:
                ticks= np.arange(si,nDays,7)
                ax[n].plot(thisCFR[si:],'-',label=place['provinces'][nn])
                ax[n].set_xticklabels(dates[ticks],{'fontsize':8})
                for label in ax[n].get_xticklabels():
                    label.set_rotation(0)
                    label.set_horizontalalignment('center')
                    label.set_fontsize(8)
        avgCFR=100*place['cfr'].mean(0)
        if move2start:
            ax[n].plot(avgCFR[si:],'k-',alpha=1, lw=3,label='%s Average CFR from provinces'%place['Country Name'])
            ax[n].plot(100*place['cfrTotCases'][si:],'k-',alpha=0.35, lw=5,label='CFR from total cases')
        else:
            ax[n].plot(avgCFR,'k-',alpha=1, lw=3,label=place['provinces'][nn])
            ax[n].plot(100*place['cfrTotCases'],'k-',alpha=0.35, lw=5,label='CFR from total cases')
        ax[n].set_ylim(0,10);
        ax[n].legend(ncol=6,loc='upper center',fontsize=8)
    figu.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.95,wspace=0.2,hspace=0.25)
    gr.ion(); gr.draw(); gr.show()
    return figu

places = [China,Australia,UK]
#
figu=plotCFRTS_Provinces(places,dates,move2start=0)
strCFR='tsam_COVID19_cfr_ProvincesChinaUKAustralia_fromFirstReport.png'
figu.savefig('./'+strCFR)
# ---------------------------------------------
strProvinces="""The

"""
