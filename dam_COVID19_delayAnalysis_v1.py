import sys
sys.path.insert(1, './')
from dam_COVID19_loadDataCSSEGIS import *
import plotly.graph_objects as go
import plotly.express as px
from ipywidgets import widgets
from bokeh.plotting import figure, output_file, show

# -------------------
# Describe the order of appearance of cases by day
# -------------------
dFC=findFirstCaseDates(totCases)
dFD=findFirstCaseDates(totDeathCases)
dFR=findFirstCaseDates(totRecovCases)
iArrival = dFC.argsort()
for n in range(37):
    print(countries[n],dFC[n],dFD[n],dFR[n],'rank: %dth'%iArrival[n])
#
d0SortedCountries = countries[iArrival]
sdFC = dFC[iArrival]
sdFD = dFD[iArrival]
sdFR = dFR[iArrival]
#
def countReportingCountriesByDay(dFC,dFD,dFR):
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
    return  nCountriesFCByDate,nCountriesFDByDate,nCountriesFRByDate

nCountriesFCByDate,nCountriesFDByDate,nCountriesFRByDate = countReportingCountriesByDay(dFC,dFD,dFR)
cumC=nCountriesFCByDate.cumsum()
cumD=nCountriesFDByDate.cumsum()
cumR=nCountriesFRByDate.cumsum()

# Delays from cumulative counts of countries reporting first cases, recoveries, and deaths
v1 = 40; v2=150
delaysCumFDC=findHDistance(days,cumC,cumD,v1,v2,nPts=v2-v1)
delaysCumFRC=findHDistance(days,cumC,cumR,v1,v2,nPts=v2-v1)
print(delaysCumFDC);  print(delaysCumFRC)
mm = findMode(delaysCumFDC)
MM = findMode(delaysCumFRC)
delayedD = cumD[mm:]
delayedR = cumR[MM:]
adjustedCD = cumC[:-mm]
adjustedCR = cumC[:-MM]
print(len(delayedD),len(adjustedCD), len(delayedR),len(adjustedCR))
#
qs = np.arange(0,1,0.01)
delaysFDC = sdFD-sdFC
delaysFRC = sdFR-sdFC
qDelaysFDC = np.quantile(delaysFDC,qs)
print('Quantiles for delays D-C:',qDelaysFDC)
medianFDCdelay = np.int32(np.median(delaysFDC) )
print('The median delay between first death and first case is %d days'%medianFDCdelay)
medianFRCdelay = np.int32(np.median(delaysFRC) )
print('The median delay between first recovery and first case is %d days'%medianFRCdelay)
delayFDFC=findMode(dFD-dFC)
delayFRFC=findMode(dFR-dFC)
print('The mode of the distribution for first death and case delays is %d'% delayFDFC)
print('The mode of the distribution for first recovery and case delays is %d'% delayFRFC)
print(np.sort(sdFD-sdFC))
print(np.sort(sdFR-sdFC))
#
# -----------------------------------
"""Poincaré diagram to plot delays in reports from cumulative counts
of first case to first recovery or death
"""
# -----------------------------------
# prepare some data
delay=medianFDCdelay
iStart = np.where(cumC>v1)[0].min()
iStop = np.where(cumC>v2)[0].min()
x=cumC[iStart:iStop]
y=cumD[iStart:iStop]
dd0= days[iStart:iStop]
#
gDelay= gr.figure(figsize=(7,5)); gr.ioff()
title="Delay between first case and first death"
gDelay.suptitle(title)
ax = gDelay.add_subplot(111)
ax.plot(days, cumC, 'b', lw=1, alpha=1, label=r'First case report')
ax.plot(days+medianFDCdelay-5, cumC, 'k', lw=2, alpha=0.4)
ax.plot(days+medianFDCdelay+5, cumC, 'k', lw=2, alpha=0.4)
ax.plot(days, cumD, 'k', lw=1, alpha=1, label=r'First death report')
ax.legend()
gr.ion(); gr.draw();gr.show()
fName="./figures_COVID19_dataAnalysis/delaysCumFDC.png"
gDelay.savefig(fName)
