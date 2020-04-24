import sys
sys.path.insert(0, './')
from NCW_baseCode import *
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


p = {'N0':120**6, 'I0':np.array([0,1,2,0]),'W0':0,'transmissionGExposure':0.5,
'exposure_N':0.455, 'exposure_I':np.array([1,1,0.1,0]),
'weight_I':np.array([0.3,0.55,0.1,0.05]), 'nI':4,
'stagesI':['Non severe','Mild/Moderate','Severe','Fatal'],
'stagesISpanish':['No graves','Moderados','Severos','Fatales'],
'infectiousTime':[14,21,25,25],'timeStep':1.0,'nSteps':360,
'pDeath':np.array([0.00001,0.0001,0.0005,0.001]),
'offset':-40, 'underReport':12}

p['timeStamps']= np.arange(0,10,)
p['tGe'] = p['timeStep']* p['transmissionGExposure']
p['pWithdraw'] = p['timeStep'] / np.array(p['infectiousTime'])
p['iC'] = np.array([p['N0'],p['W0'],p['I0']])
susc= np.int32( p['N0'] * p['exposure_N'] * p['weight_I'])
newI= binom.rvs( susc[0], p['tGe'] * np.dot( p['exposure_I'], p['I0'] ),loc=0, size=1)
randNIWIncidence(Z=p['iC'], p=p)

MexicoCases_April16 =np.array ([3, 4, 5, 6, 7, 11, 15, 26, 41, 53, 82, 93, 118, 164, 203, 251, 316, 367, 405, 475, 585, 717, 848, 993, 1094, 1215, 1378, 1510, 1688, 1890, 2143, 2439, 2785, 3181, 3441, 3844, 4219, 4661, 5014, 5399, 5847, 6297])
MexicoDeaths_April16=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 2, 4, 5, 6, 8, 12, 16, 20, 28, 29, 37, 50, 60, 79, 94, 125, 141, 174, 194, 233, 273, 296, 332, 406, 449, 486])

N,W,I= randNIWDynamics(p, rhs=randNIWIncidence, nonAutPars={})


def plotEpidemic(N,W,I,p,maxCases=10000,maxDeaths=1000):
    INon= np.array([I[m][0] for m in range(len(I))])
    IMild= np.array([I[m][1] for m in range(len(I))])
    ISevere= np.array([I[m][2] for m in range(len(I))])
    IFatal= np.array([I[m][3] for m in range(len(I))])
    ITot = INon+IMild+ISevere+IFatal
    dead = list()
    for nn in range(p['nI']):
        dead.append(p['pDeath'][nn]*W)
    days = np.arange(0,len(N)) +p['offset']
    daysMexico = np.arange(0,len(MexicoCases_April16))
    peakInd =ITot.argmax()
    peakDay = days[peakInd]

    fig= gr.figure(figsize=(11,5)); gr.ioff()
    rows=1;cols=1; ax=list()

    for rc in range(rows*cols):
        ax.append(fig.add_subplot(rows,cols,rc+1))

    cax=inset_axes(parent_axes=ax[0],
                            width="30%", # width = 30% of parent_bbox
                            height="40%", # height : 1 inch
                            loc='center right')
    #dax=inset_axes(parent_axes=ax[2],
    #                        width="30%", # width = 30% of parent_bbox
    #                        height="50%", # height : 1 inch
    #                        loc='center left')
    #ax[1].plot(days, N,label='$N$ (no infectados)')
    #ax[1].plot(days, W,label='$W$ (recuperados)')
    ax[0].plot(days, INon,label=r'$I_{Non}$')
    ax[0].plot(days, IMild,label=r'$I_{Mild}$')
    ax[0].plot(days, ISevere,label=r'$I_{Severe}$')
    ax[0].plot(days, IFatal,label=r'$I_{Fatal}$')
    ax[0].plot(days, ITot,label=r'$Tot Cases$')
    str1=r'Confirmed cases, S. Salud Mexico (04/16/2020)'
    ax[0].plot(daysMexico, p['underReport']*MexicoCases_April16,'o',ms=2,label=str1)
    cax.plot(days, INon)
    cax.plot(days, IMild)
    cax.plot(days, ISevere)
    cax.plot(days, IFatal)
    cax.plot(days, ITot)
    cax.plot(daysMexico, p['underReport']* MexicoCases_April16,'o',ms=2,)
    #for rc in range(p['nI']):
    #    ax[2].plot(days, dead[rc],label=r'$Fallecidos_{%s}$'%p['stagesISpanish'][rc])
    #ax[2].plot(days, dead[rc].sum(0),label=r'$Fallecidos_{%s}$'%p['stagesISpanish'][rc])
    #dax.plot(days, dead[rc].sum(0))
    #ax[2].plot(daysMexico, underReport*MexicoDeaths_April16,label=r'$Fallecidos_{04/16}$')
    #dax.plot(daysMexico, underReport*MexicoDeaths_April16,label=r'$Fallecidos_{04/16}$')
    for rc in range(rows*cols):
        ax[rc].legend()
    cax.set_xlim(0,50); cax.set_ylim(0,maxCases)
    ax[0].set_ylabel('Confirmed cases')
    ax[0].set_xlabel('Days from first report on Feb 27, 2020')
    #dax.set_xlim(0,50); dax.set_ylim(0,maxDeaths)
    strTit= '''Covid-19 dynamics from the contribution of infectious individuals with different contagion intervals \n
    (Herrera-Nolasco, Herrera-Valdez, 2020)'''
    ax[0].text(peakInd+5, ITot.max()*.95, '%d days, %d infected'%(peakDay,ITot[peakInd]))
    ax[0].text(peakInd+5, ITot.max()*.9, 'Peak estimated around May 13, 2020')
    fig.subplots_adjust(left=0.075,bottom=0.075,right=0.97,top=0.85,wspace=0.2,hspace=0.25)
    fig.suptitle(strTit)
    gr.ion(); gr.draw(); gr.show()
    return fig

p['offset']=-40
fig = plotEpidemic(N,W,I,p,maxCases=40000,maxDeaths=2000)
