# tsam_COVID19_cfr_JHU
from tsam_COVID19_baseCode import *
cwd = os.getcwd()


cases,deathCases,recovCases = getCSSEGISandData(urlData=1)
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
gdeathCases,countries_deathCases = sortDataByCountry(deathCases, nHeaderCols)
gRecov,countries_Recov = sortDataByCountry(recovCases, nHeaderCols)

# -------------------
# Sum the counts from each country and construct a new array
# -------------------
# These arrays have the same size as the countries array (unique countries)
wcCases=gatherDataByCountry(df=cases,nHeaderCols=4)
wcDeathCases=gatherDataByCountry(df=deathCases,nHeaderCols=4)
wcRecovCases=gatherDataByCountry(df=recovCases,nHeaderCols=4)

# -------------------
# Delays among the first cases
# -------------------

# Order of appearance of cases from the reports
iFC=findFirstCaseDates(wcCases)
iFD=findFirstCaseDates(wcDeathCases)
iFR=findFirstCaseDates(wcRecovCases)
iArrivalReports = iFC.argsort()
d0SortedCountries = countries[iArrivalReports]
siC = iFC[iArrivalReports]
siD = iFD[iArrivalReports]
siR = iFR[iArrivalReports]
print('%First case, & to first death & to first recovery \\\\')
print('\hline &&&\\\\')
for mm in range(nCountries):
    #print('%s & %s & %s & %s\\\\'%(countries[jj[mm]], dates[siC[mm]], dates[siD[mm]], dates[siR[mm]]))
    #print('& %d & %d & %d \\\\ '%(siC[mm],siD[mm]-siC[mm],siR[mm]-siC[mm]))
    sss = [d0SortedCountries[mm], dates[siC[mm]], siC[mm], dates[siD[mm]], siD[mm]-siC[mm], dates[siR[mm]], siR[mm]-siC[mm]]
    print('{s[0]} & {s[1]} ({s[2]})  & {s[3]} ({s[4]})  & {s[5]} ({s[6]}) \\\\'.format(s=sss))
    if (mm<nCountries-1):
        if (siC[mm]<siC[mm+1]):
            print('\hline ')

# -------------------------
# Time Course of the reports
# Quantify the arrival of the reports
# -------------------------
casesCountriesByDate= list()
deathsCountriesByDate= list()
recovsCountriesByDate= list()
for m in range(nDays):
    casesCountriesByDate.append(len(np.where(siC==m)[0]))
    deathsCountriesByDate.append(len(np.where(siD==m)[0]))
    recovsCountriesByDate.append(len(np.where(siR==m)[0]))
casesCountriesByDate = np.array(casesCountriesByDate)
deathsCountriesByDate = np.array(deathsCountriesByDate)
recovsCountriesByDate = np.array(recovsCountriesByDate)

#
ticks= np.arange(0,nDays,7)
fRepTS= gr.figure(figsize=(11,5)); rows=1; cols=2
ax=fRepTS.add_subplot(rows,cols,1);
bx=fRepTS.add_subplot(rows,cols,2);
gr.ioff()
bx.scatter(siC,siD-siC,s=10, c='r', marker='o')
bx.scatter(siC,siR-siC,s=10, c='g', marker='o')
bx.set_ylabel('Days from first case reported')
bx.set_xlabel('Days after $d_0$')
ax.plot(dates, casesCountriesByDate.cumsum(),label='Cases')
ax.plot(dates, deathsCountriesByDate.cumsum(),label='Deaths')
ax.plot(dates, recovsCountriesByDate.cumsum(),label='Recovs')
ax.set_xticks(ticks)
ax.set_xticklabels(dates[ticks])
ax.legend(loc='upper left')
for label in ax.get_xticklabels():
    label.set_rotation(90)
    label.set_horizontalalignment('center')
    label.set_fontsize(8)
gr.ion(); gr.draw()

n=10; print(r'%d countries reported within the first %d days from $d_0$'%(len(siC[siC<n+1]),n))
countries[jj[siC<n+1]]

"""
The first countries that reported cases of COVID-19 on January 22, 2020 ($d_0$)
were South Korea, China, Taiwan, US, Japan, and Thailand.
Australia reported their first case 4 days later, together with Canada, on January 26, 2020.
The first reports from Europe arrived 5 days later from Germany, on January 27.
The first report in the Middle East came from the United Arab Emirates on January 29, 7 days after the first case reported.
Therefore, the COVID-19 epidemic had been reported from all continents but Africa one week after the first case was reported, and the first report from Africa arrived 23 days later from Egypt, on February 14, 2020. These data indicate that the COVID-19 epidemic was wide spread throughout the globe within 2 weeks of the first case, strongly suggesting that the dissemination of the SARS-CoV-2 virus among the human population occurred with a delay of approximately 2 weeks.

These data reflect the delay in reports from different countries relative to 𝑑0
. Of note, 25 countries reported cases within 10 days from 𝑑0, and 31 countries had reports within the first 30 days. Of interest to the authors, Mexico reported their first case 37 days later, on February 28, 2020, together with Nigeria, New Zealand, Lithuania, Iceland, and Belarous. The last country that reported cases was Sao Tome and Principe, on April 6th, 2020, 75 days after the first report.
"""
# --------------------------------------------
# Histograms of the first arrivals for cases, deaths, recoveries
# --------------------------------------------
cCases, bins = np.histogram(iFC,np.arange(0,iFC.max()))
cDeathCases, bins = np.histogram(iFD-iFC,np.arange(0,iFC.max()))
cRecov, bins = np.histogram(iFR-iFC,np.arange(0,iFC.max()))
#
cCases, binsC = np.histogram(iFC,np.arange(0,(iFD-iFC).max()))
cDeathCases, binsD = np.histogram(iFD-iFC,np.arange(0,(iFD-iFC).max()))
cRecovCases, binsR = np.histogram(iFR-iFC,np.arange(0,(iFD-iFC).max()))
#
def plotReportingDelays():
    fDelays= gr.figure(figsize=(7,10)); gr.ioff(); rows=3; cols=1;
    ax=list(); tAx=list()
    fDelays.suptitle('Delays in first case report, first death, and first recovery, all countries')
    for m in range(rows*cols):
        ax.append(fDelays.add_subplot(rows,cols,m+1))
        tAx.append(ax[m].twinx())
    tAx[0].plot(casesCountriesByDate.cumsum(),'b',lw=2, label='Cumulative # countries reporting cases')
    tAx[0].plot(deathsCountriesByDate.cumsum(),'k',lw=2,label='Cumulative # countries reporting deaths')
    tAx[0].plot(recovsCountriesByDate.cumsum(),'orange',lw=2, label='Cumulative # countries reporting recoveries')
    ax[0].bar(binsC[:-1],cCases,width=0.8)
    ax[1].bar(binsD[:-1],cDeathCases,width=0.8,label='Death')
    ax[2].bar(binsR[:-1],cRecovCases,width=0.8,label='Recoveries')
    ax[0].set_xlabel('Days from $d_0$')
    ax[1].set_xlabel('Days from first case reported (same country)')
    ax[2].set_xlabel('Days from first case reported (same country)')
    for m in range(rows*cols):
        ax[m].set_yticks(np.arange(0,np.maximum(cRecovCases.max()+1,cCases.max()),2))
        ax[m].set_ylabel('# Countries',rotation=90)
        tAx[m].legend(loc='upper left')
    fDelays.subplots_adjust(left=0.075,bottom=0.075,right=0.9,top=0.95,wspace=0.1,hspace=0.25)
    gr.ion(); gr.draw(); gr.show()
    fDelays.savefig('tsam_COVID19_JHU_delaysAllCountries.png')
    return fDelays

fDelays= plotReportingDelays()
n=7; print(r'%d countries reported within the first %d days from $d_0$'%(len(siC[siC<n+1]),n))
countries[jj[(siD-siC) <n+1]]
strDelays = '''
The delays between the first case reported and the first death reported within a single country have a wide distribution, between 0 and 76 days.
However, the bulk of the distribution can be found between 0 and 30 days, with a mode at 16 days. In consideration of the time course of the infection by SARS-CoV-2 reported in the literature \citep{}, the delay between death and contagion can be 21 or more days \citep{}.

Assuming that the testing started at least one week before the date of the first case reported, the data suggest that the epidemic in those countries contributing to left portion of the histogram had started at least 3 weeks before their report. The mode at 16 days fits within the 5-7 days incubation periods reported plus the 2 weeks to develop severe symptoms. For those countries that reported their first death at the mode or after, these data suggest that the epidemic started in those countries at least one week before their first case reported.

The delays between the first recovery and the first case reported within a single country are distributed between 0 and 38 days. However, for most countries this delay is between 0 and 27 days approximately, with two clear modes at 12 and 16 days.
'''
# -------------------
# Search regions to illustrate the case-fatality ratios
# -------------------
Pops_Millions = {'China':1439323776, 'Japan':126476461,'Korea, South':51269185, 'Indonesia':273523615, 'Singapore':5850342, 'Mexico':128932753,
'US':331002651, 'Canada':37742154, 'Argentina':45195774, 'Brazil':212559417, 'Colombia':50882891, 'Niger':24206644, 'Algeria':43851044, 'Egypt':102334404, 'South Africa':59308690,'Spain':46754778,'Italy':60461826,'France':67886011,'Germany':83783942, 'Australia':25499884,'United Kingdom':67886011,'Iran':83992949,'Israel':8655535}
#
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

#
def plotCasesDeathsTS(casesTS,deathsTS,regions,countries, convFactor=1000,saveFig=1):
    ii = getIndsRegions(countries, regions)
    figu= gr.figure(figsize=(7,9))
    if convFactor <= 1000:
        figu.suptitle('Deaths vs cases per {: d} habitants'.format(convFactor))
    elif convFactor == 10**6:
        figu.suptitle('Deaths vs cases per million ')
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
            cas=convFactor*np.float64(casesTS[region[nn]])/Pops_Millions[countries[region[nn]]]
            dea=convFactor*np.float64(deathsTS[region[nn]])/Pops_Millions[countries[region[nn]]]
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
    if saveFig>0:
        figName='tsam_COVID19_JHU_cases-deaths_x%d_JHU.png'%convFactor
        figu.savefig('./'+figName)
    return figu

#
ii = getIndsRegions(countries, regions)
figu=plotCasesDeathsTS(casesTS=totCases,deathsTS=totDeathCases,regions=regions,countries=countries, convFactor=10**6)




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
cfr= correctedArrayRatio(wcDeathCases,wcCases)
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
China['cfrwcCases']=cfr[np.where(countries=='China')[0][0]]
China['cases'], China['indsCases']= gatherDataSingleCountry(cases,'China')
China['deathCases'], China['indsdeathCases']= gatherDataSingleCountry(deathCases,'China')
China['cfr'] = correctedArrayRatio(China['deathCases'],China['cases'])
China['provinces']=cases.iloc[China['indsCases'],0].to_numpy()
China['nProvinces'] = len(China['provinces'])
China['startDaysCases']= findCaseStarts(places=China['provinces'],cases=China['cases'])
#
UK =  {'Country name':'UK'}
UK['cfrwcCases']=cfr[np.where(countries=='United Kingdom')[0][0]]
UK['cases'], UK['indsCases']= gatherDataSingleCountry(cases,'United Kingdom')
UK['deathCases'], UK['indsdeathCases']= gatherDataSingleCountry(deathCases,'United Kingdom')
UK['cfr'] = correctedArrayRatio(UK['deathCases'],UK['cases'])
UK['provinces']=cases.iloc[UK['indsCases'],0].to_numpy()
UK['nProvinces'] = len(UK['provinces'])
UK['startDaysCases']= findCaseStarts(places=UK['provinces'],cases=UK['cases'])
UK['provinces'][6]= 'Great Britain'
print(UK['provinces'])
#
Australia =  {'Country name':'Australia'}
Australia['cfrwcCases']=cfr[np.where(countries=='Australia')[0][0]]
Australia['cases'], Australia['indsCases']= gatherDataSingleCountry(cases,'Australia')
Australia['deathCases'], Australia['indsdeathCases']= gatherDataSingleCountry(deathCases,'Australia')
Australia['cfr'] = correctedArrayRatio(Australia['deathCases'],Australia['cases'])
Australia['provinces']=cases.iloc[Australia['indsCases'],0].to_numpy()
Australia['nProvinces'] = len(Australia['provinces'])
Australia['startDaysCases']= findCaseStarts(places=Australia['provinces'],cases=Australia['cases'])
#
def plotCFRTS_Provinces(place,dates,move2start=1):
    figu = gr.figure(figsize=(13,7))
    figu.suptitle('Percentage of dead/confirmed between %s-%s in %s'''%(dates[0],dates[-1],place['Country name']))
    gr.ioff();cols=1; rows = 1
    ax=figu.add_subplot(rows,cols,1)
    print(place['Country name'])
    ticks= np.arange(0,nDays,7)
    si = place['startDaysCases'].min()
    for nn in range(len(place['provinces'])):
        #print(place['provinces'][nn])
        thisCFR=100*place['cfr'][nn]
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
    avgCFR=100*place['cfr'].mean(0)
    if move2start==1:
        ax.plot(avgCFR[si:],'k-',alpha=1, lw=3,label='Average CFR from provinces in %s'%place['Country name'])
        ax.plot(100*place['cfrwcCases'][si:],'k-',alpha=0.35, lw=5,label='CFR from total cases')
        strCFR='tsam_COVID19_JHU_cfr_Provinces'+place['Country name']+'_fromDate0.png'
    else:
        ax.plot(avgCFR,'k-',alpha=1, lw=3,label='Average CFR from provinces in %s'%place['Country name'])
        ax.plot(100*place['cfrwcCases'],'k-',alpha=0.35, lw=5,label='CFR from total cases')
        strCFR='tsam_COVID19_JHU_cfr_Provinces'+place['Country name']+'_fromFirstLocalReport.png'
    ax.set_ylim(0,13);
    ax.legend(ncol=6,loc='upper center',fontsize=8)
    figu.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.95,wspace=0.2,hspace=0.25)
    figu.savefig('./'+strCFR)
    gr.ion(); gr.draw(); gr.show()
    return figu
#
figuChina=plotCFRTS_Provinces(China,dates,move2start=0)
figuUK=plotCFRTS_Provinces(UK,dates,move2start=0)
figuAust=plotCFRTS_Provinces(Australia,dates,move2start=0)
# ---------------------------------------------

# ----------------------------------
# Age cfrs and inference about probability of death due to COVID-19
# ----------------------------------
SKorea={'ageCFR': np.array([0,0,0,0.001,0.001,0.004,0.015,0.063,0.116]),'CFR_2020?':0.02,'Country name':'Korea, South'}
Italy={'ageCFR': np.array([0,0,0,0.003,0.004,0.01,0.035,0.128,0.202]),'CFR_20200317':0.072,'Country name':'Italy'}
China['ageCFR'] = np.array([0,0.002,0.002,0.002,0.004,0.013,0.036,0.08,0.14])
China['CFR_20200211'] = 0.02
#
aGroups=np.arange(0,90,10)
ageCountries = [China,Italy,SKorea]
ages = np.arange(0,90)
#
f11=gr.figure(figsize=(11,5)); gr.ioff(); rows =1;cols=1
f11.suptitle('CFR by age')
lshift=2.5;  ax=f11.add_subplot(rows,cols,1)
for n in range(len(ageCountries)):
    ax.bar(aGroups+ (n*lshift), 100*ageCountries[n]['ageCFR'],width=2,align='edge',label=ageCountries[n]['Country name']);
ax.legend(loc='upper left')
ax.set_xticks(aGroups)
ax.set_yticks(np.arange(0,21,2))
ax.set_xlabel('Age groups')
gr.ion(); gr.draw();  gr.show()
f11Name='tsam_COVID19_JHU_cfr+propDeathCasesByAge_China+SKorea+Italy.png'
f11.savefig(f11Name)


# -------------------------------------------
# Calculation for proportions of death by age group for different ages and prop deathCases by age group using Bayes' theorem
# -------------------------------------------
strCPropsDead='''
We can establish qualitative estimations of the percentage of deaths based on these age CFRs. To do so, we calculate the weights of the CFRs in each age group, which represent the proportion of deaths in each group.

#### Calculation of proportions of death within cases, by age groups

Let $c_i$ represent the case-fatality ratio in $i$th age group, $i \in \left\{0,1,...,8\right\}$
The weight of the $i$th age group on the total death reported via CFRs is
$$w_i = \frac{c_i}{\sum_{j=0}^{n-1}} w_j.$$
Then the probability of death in the $i$th age group can be estimated by multiplying the CFR and the weight of the $i$th age group on the CFR. Then, multiplication by the number of cases would yield the death toll by age group.

A few sutile differences can be noted from these data. For instance, the CFRs for the first 3 groups in Italy and South Korea are very similar. In Italy, the percentage of dead cases is very large among the last two age groups, which in part is due to the fact that Italians implemented a triage system that sometimes denied respirators to patients with small chances of survival.
'''

for n in range(len(ageCountries)):
    C = ageCountries[n]
    print(C['Country name'])
    ageCFRweights = C['ageCFR'].sum()
    C['ageDeadCaseProps'] = C['ageCFR']/ageCFRweights
    C['cfrwcCases']=cfr[np.where(countries==C['Country name'])[0][0]]
    C['wcCases'] = wcCases[np.where(countries==C['Country name'])[0][0]]
    C['wcDeathCases'] = wcDeathCases[np.where(countries==C['Country name'])[0][0]]
    C['wcDeathCases_ageProps'] = np.zeros((len(aGroups),len(C['wcDeathCases'])),'int64')
    for nn in range(len(aGroups)):
        C['wcDeathCases_ageProps'][nn,:]= C['ageDeadCaseProps'][nn] * C['wcCases']* C['cfrwcCases']
        print(C['wcDeathCases_ageProps'][nn])


a0s = [74,70,72]; aMaxs=[1.4*0.12,1.3*0.202,1.3*0.14]
pMaxs=[0.76,0.68,0.68]
f10=gr.figure(figsize=(13,11)); gr.ioff();
f10.suptitle('Comparison between CFRs by age and proportion of deaths by age (among confirmed)')
axCFR=list(); axDP=list(); rows =len(ageCountries);cols=2
for n in range(len(ageCountries)):
    C=ageCountries[n]
    axCFR.append(f10.add_subplot(rows,cols,2*n+1))
    axCFR[n].text(0,0.1,C['Country name'])
    axCFR[n].bar(aGroups, C['ageCFR'],width=8,align='edge',label='CFR '+C['Country name']);
    axCFR[n].plot(ages, sigmoid(a=ages,a0=a0s[n],aMax=aMaxs[n],n=9))
    axCFR[n].legend(loc='upper left')
    axDP.append(f10.add_subplot(rows,cols,2*n+2))
    axDP[n].text(0,0.1,C['Country name'])
    axDP[n].bar(aGroups, C['ageDeadCaseProps'],width=8,align='edge',label='P(COVID-19 dead by age) '+ageCountries[n]['Country name']);
    axDP[n].plot(ages, sigmoid(a=ages,a0=a0s[n],aMax=pMaxs[n],n=9))
    axDP[n].legend(loc='upper left')
gr.ion(); gr.draw();gr.show()
#
f10Name='tsam_COVID19_JHU_cfr+propDeathCasesByAge_China+SKorea+Italy.png'
f10.savefig(f10Name)



# ------------------------------------------
# Comparison of death estimates by age and deaths reported
# ------------------------------------------
ages = np.arange(0,90)
f12=gr.figure(figsize=(13,9)); gr.ioff(); rows =3;cols=1
f12.suptitle('P(death) due to COVID-19 by age')
ax=list()
ticks= np.arange(0,nDays,7)
convFactor =1
for n in range(len(ageCountries)):
    ax.append(f12.add_subplot(rows,cols,n+1))
    C = ageCountries[n]
    print(C['Country name'])
    ax[n].plot(C['wcDeathCases'],'k--', alpha=1, lw=1,label='Deaths reported '+ C['Country name']);
    ax[n].plot(C['wcDeathCases_ageProps'].sum(0),'k',lw=1,label='Sum of estimates '+ C['Country name']);
    for nn in range(len(aGroups)):
        ax[n].plot(C['wcDeathCases_ageProps'][nn],label=C['Country name']+ '[%d,%d)]'%(aGroups[nn],aGroups[nn]+10));
    ax[n].legend(loc='upper left',ncol=3, fontsize=8)
    ax[n].set_xticks(ticks)
    ax[n].set_xticklabels(dates[ticks],{'fontsize':8})
    for label in ax[n].get_xticklabels():
        label.set_rotation(0)
        label.set_horizontalalignment('center')
        label.set_fontsize(8)
gr.ion(); gr.draw();  gr.show()

f12Name='tsam_COVID19_JHU_cfr+propDeathCasesByAgeTS_China_SKorea_Italy.png'
f12.savefig(f12Name)

# -----------------------------------------------
# Estimation of cases by age in Mexico
# -----------------------------------------------
Mexico = dict()
Mexico['cfr'] =cfr[np.where(countries=='Mexico')[0][0]]
Mexico['wcCases']=wcCases[np.where(countries=='Mexico')[0][0]]
Mexico['wcDeathCases']=wcDeathCases[np.where(countries=='Mexico')[0][0]]
Mexico['startInd'] = np.where(Mexico['wcCases']>0)[0].min()
Mexico['popMillions']= Pops_Millions['Mexico']
Mexico['delayFromReport0']=np.where(Mexico['wcCases']>0)[0].min()
print('First case reported %d days after Report 0'%Mexico['delayFromReport0'])

convFactors=[1,8,12,16];

def estimateDeathsByAgeMexico(Mexico, dates, subReportFactor=1):
    si = Mexico['delayFromReport0']
    ticks= np.arange(0,len(dates),7)
    casesMexico = subReportFactor * Mexico['wcCases']
    #pConfirmation = casesMexico/Mexico['popMillions']
    f13=gr.figure(figsize=(13,9)); gr.ioff(); rows =3;cols=1;ax=list()
    for n in range(len(ageCountries)):
        ax.append(f13.add_subplot(rows,cols,n+1))
        C = ageCountries[n]
        ax[n].plot(subReportFactor* Mexico['wcDeathCases'],'k',lw=1,label='Deaths reported in Mexico');
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
        ax[n].set_xlim(si,len(Mexico['wcDeathCases']))
    f13.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.95,wspace=0.1,hspace=0.25)
    f13Name='tsam_COVID19_JHU_cfr+propDeathCasesByAgeTS_EstimatesMexico_subReportFactor%d.png'%subReportFactor
    gr.ion(); gr.draw();  gr.show()
    f13.savefig(f13Name)
    return f13

f13_1 = estimateDeathsByAgeMexico(Mexico, dates,subReportFactor=1)
f13_10 = estimateDeathsByAgeMexico(Mexico, dates,subReportFactor=10)
f13_12 = estimateDeathsByAgeMexico(Mexico, dates,subReportFactor=12)
