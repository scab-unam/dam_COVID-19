import pylab as gr
import matplotlib.pyplot as plt
import scipy as sc
sc.test('all')


def sir(U,t,beta,mu):
    S,I,R=U
    N=S+I+R
    newInfections= beta*S*I/N
    newRecoveries= mu*I
    dS = -newInfections
    dI = newInfections-newRecoveries
    dR = newRecoveries
    return dS, dI, dR

def uvw(Z,t,beta,delta):
    u,v,w=Z
    du = -beta*u*v
    dw = delta*v
    dv = -du -dw
    return du, dv, dw

def uvField(ax,beta=0.5, delta=0.1, parts=100j):
    V,U=sc.mgrid[0:1:parts, 0:1:parts]
    dw=delta*V
    du= -beta*U*V
    dv= -du - dw
    #speed = sc.sqrt( U*U + V*V)
    speed = 2*sc.sqrt(V*V)
    #lw=1; 
    db=delta/beta
    print('delta / beta = %f'%(db))
    lw=2*speed/speed.max()
    ax.plot([db,db],[0,1],'k--',lw=1)
    #ax.streamplot(U,V, du,dv, color='k', linewidth=lw, cmap=gr.cm.gray_r)
    #ax.streamplot(U,V, u,v, color=v, linewidth=2, cmap=gr.cm.autumn)
    #ax.streamplot(U,V, u,v, color='k', density=[0.75, 0.75], linewidth=lw)
    ax.streamplot(U,V, du,dv, color='k', linewidth=lw)
    return V,U,du,dv


anhos= sc.arange(1990,2012.1,1)
casosSIDA= sc.array([10185,13548,17383,21296,26009,31062,36638,42933,49799,58565,67251,75754,84081,92296,100596,109308,117783,125496,132712,139722,147045,152685,158302])

def vLogistic(v,t,beta,delta):
    dv = beta * v * (1 - delta/beta - v)
    return dv

def vLogisticBack(v,t,beta,delta):
    dv = - beta * v * (1 - delta/beta -v)
    return dv


#def ajusteLogisticoVIH(ax, 
vinfs=sc.arange(0.0015,0.0025,0.0001)#, 
betas=sc.arange(1,130,1.0)#, v0=1e-6, 
popSize=110e6#, tmax=50.0, tstep=0.001, 
startYear=1960.0#): 
casos=casosSIDA/popSize
if 0:
    vStart= casos[0]
    beta=betas[0]
    delta=(1-vinfs[0])*beta
    sampTime=sc.arange(anhos[0], anhos[-1]+1, 0.001)
    residue=1e13
    for b in betas:
        for v in vinfs:
            d=(1-v)*b
            #orbit=sc.integrate.odeint(vLogisticBack,casos[-1],sampTime,args=(b,d),rtol=1e-6).transpose()
            orbit=sc.integrate.odeint(vLogistic,casos[0],sampTime,args=(b,d),rtol=1e-6).transpose()
            #vOrbit=sc.hstack([orbit[0][-1:0:-1], orbit[0][0]])
            vOrbit=orbit[0]
            thisResidue=0.0
            for m in range(0,len(anhos)):
                thisIndex=gr.find(sampTime>anhos[m]).min() -1
                thisResidue= thisResidue+ (vOrbit[thisIndex] - casos[m])**2
            
            if residue>thisResidue:
                print('(%g, %g)'%(residue,thisResidue)) 
                residue=thisResidue
                beta = b; delta=d
                print('(new b, new d, new vinf)=(%g, %g, %g)'%(beta,delta, v))
    
    vinf = 1- delta/beta
    print('(beta, delta, vinf)=(%g, %g, %g)'%(beta,delta, vinf))
    gr.ion()
    gr.figure()
    gr.plot(sampTime, vOrbit, 'k', lw=1, alpha=1.0)
    gr.plot(anhos, casos, 'wo', ms=3, alpha=1.0)
    sampTime=sc.arange(anhos[0], 2020.0, 0.001)
    orbitF=sc.integrate.odeint(vLogistic,vStart,sampTime,args=(beta,delta),rtol=1e-6).transpose()
    sampTime=sc.arange(anhos[0], anhos[0]-startYear, 0.001)
    orbitB=sc.integrate.odeint(vLogisticBack,vStart,sampTime,args=(beta,delta),rtol=1e-6).transpose()
    vOrbit=sc.hstack([orbitB[0],orbitF[0]])
    sampTime=sc.arange(anhos[0], 2020.0, 0.001)
    if 0: 
        ax.plot(sampTime, vOrbit, 'k', lw=1, alpha=1.0)
        ax.plot(anhos, casos, 'wo', ms=3, alpha=1.0)
        ax.set_ylabel(u'% casos')
        ax.set_xlabel(u"años")
    #return vOrbit



def seasonalHumps(t,centers,amp=1.0):
    ind= amp*sc.exp( -(t-centers)**2 / 2.0 ).sum()
    return ind

def sirPeriodicBeta(U,t,beta,mu,amp):
    S,I,R=U
    N=S+I+R
    newInfections= beta*betaFun(t,amp)*S*I/N
    newRecoveries= mu*I
    dS = -newInfections
    dI = newInfections-newRecoveries
    dR = newRecoveries
    return dS, dI, dR

