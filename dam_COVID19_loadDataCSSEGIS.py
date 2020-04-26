import os
import numpy as np
import pandas as pd
import datetime as dt
import sys
sys.path.insert(1, './')
from dam_COVID19_baseCode import *
import matplotlib.pylab as gr
small={'family' : 'normal','weight' : 'normal','size'   : 8}
medium={'family' : 'normal','weight' : 'normal','size'   : 10}
large={'family' : 'normal','weight' : 'bold','size'   : 13}
gr.rc('font', size=small['size'], weight='normal')          # controls default text sizes
gr.rc('axes', titlesize=small['size'])     # fontsize of the axes title
gr.rc('axes', labelsize=medium['size'])    # fontsize of the x and y labels
gr.rc('xtick', labelsize=small['size'])    # fontsize of the tick labels
gr.rc('ytick', labelsize=small['size'])    # fontsize of the tick labels
gr.rc('legend', fontsize=small['size'])    # legend fontsize
gr.rc('figure', titlesize=large['size'])  # fontsize of the figure title
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# Example:
# inset_axes = inset_axes(parent_axes,
#                     width="30%", # width = 30% of parent_bbox
#                     height=1., # height : 1 inch
#                     loc=3)
# Indexing and gathering of data

datosAbiertosCovid19Mexico='http://187.191.75.115/gobmx/salud/datos_abiertos/datos_abiertos_covid19.zip'
datosMexico=pd.read_csv(datosAbiertosCovid19Mexico, compression='zip',encoding='latin-1')
mexReference="""Data for Mexico obtained from http://187.191.75.115/gobmx/salud/datos_abiertos/datos_abiertos_covid19.zip"""
print(mexReference)
strReference="""Data obtained from https://github.com/CSSEGISandData/COVID-19/"""
print(strReference)


cases, deathCases,recovCases = getCSSEGISandData(urlData=1)
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
npCases = cases.to_numpy()
countries = np.unique(npCases[:,1])
nCountries = len(countries)
print('Considering data from {d} countries'.format(d=nCountries))
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
    report nCountriesFCByDate,nCountriesFDByDate,nCountriesFRByDate

nCountriesFCByDate,nCountriesFDByDate,nCountriesFRByDate = countReportingCountriesByDay(dFC,dFD,dFR)
cumC=nCountriesFCByDate.cumsum()
cumD=nCountriesFDByDate.cumsum()
cumR=nCountriesFRByDate.cumsum()
#
# -------------------
# Search regions to illustrate the case-fatality ratios
# -------------------
Pops_Millions = { 'Algeria':43851044, 'Argentina':45195774, 'Australia':25499884, 'Brazil':212559417, 'Bolivia':11353142, 'Canada':37742154, 'China':1439323776, 'Colombia':50882891, 'Egypt':102334404, 'France':67886011, 'Germany':83783942, 'Indonesia':273523615, 'Iran':83992949, 'Israel':8655535, 'Italy':60461826, 'Japan':126476461, 'Korea, South':51269185, 'Mexico':128932753, 'Niger':24206644, 'Singapore':5850342, 'South Africa':59308690, 'Spain':46754778, 'United Kingdom':67886011, 'US':331002651}
