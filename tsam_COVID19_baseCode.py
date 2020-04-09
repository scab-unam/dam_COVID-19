import pandas as pd
import numpy as np
import matplotlib.pylab as gr
import datetime
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# Example:
# inset_axes = inset_axes(parent_axes,
#                     width="30%", # width = 30% of parent_bbox
#                     height=1., # height : 1 inch
#                     loc=3)
# Indexing and gathering of data
def getIndsSingleRegion(countries, region):
    i = [np.where(countries==region[nn])[0][0] for nn in range(len(region))]
    return i

def getIndsRegions(countries, regions):
    ii =list()
    for mm in range(len(regions)):
        reg = regions[mm]
        ii.append(getIndsSingleRegion(countries,reg))
    return ii
#
def sortDataByCountry(df,nHeaderCols=4):
    """
    Sort data from the same country
    Output:
    x --> Sorted data
    places --> countries sorted as in x
    """
    iCoun= df.iloc[:,1].to_numpy().argsort()
    x = df.iloc[iCoun,nHeaderCols:].to_numpy()
    places = df.iloc[iCoun,1].to_numpy()
    return x, places

def gatherDataSingleCountry(df,country, nHeaderCols=4):
    """
    Gather data from the same country
    """
    cc = df.iloc[:,1].to_numpy()
    countries = np.unique(df.iloc[:,1].to_numpy())
    nCountries = len(countries);
    x = list()
    iC= np.where(cc== country)[0]
    a =df.iloc[iC,nHeaderCols:].to_numpy()
    a.sum(0)
    return a, iC


def gatherDataByCountry(df,nHeaderCols=4):
    """
    Gather data from the same country
    """
    cc = df.iloc[:,1].to_numpy()
    countries = np.unique(df.iloc[:,1].to_numpy())
    nCountries = len(countries);
    x = list()
    for n in range(nCountries):
        iC= np.where(cc== countries[n])[0]
        a =df.iloc[iC,nHeaderCols:].to_numpy()
        x.append(a.sum(0))
    return np.array(x)

def findCaseStarts(places,cases):
    """
    Find the indices of at which the first cases are observed in each location from the list places.
    """
    startInds = list()
    for n in range(len(places)):
        ii= np.where(cases[n]>0)[0]
        if len(ii)>0:
              startInds.append(ii.min())
        else:
              startInds.append(len(cases[n])-1)
    return np.array(startInds)

# -------------------
# Ratios by array. Correction in cases there are zeros in the denominators
# -------------------
def correctedArrayRatio(a,b):
    """
    correctedArrayRatio(a,b) calculates a/b from arrays a and b,
    assuming that a is less than or equal to b,
    and correcting for possible zeros in b.
    In those cases, the value of b is set to 1, and the quotient is still 0
    """
    bCorr= b.copy()
    bCorr[bCorr==0]=1
    return np.float64(a)/np.maximum(b,bCorr)

def add_subplot_axes(ax,rect):
    fig = gr.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax
