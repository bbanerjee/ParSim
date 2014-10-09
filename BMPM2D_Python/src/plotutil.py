import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy as copy
from matplotlib.widgets import Slider, Button, RadioButtons

def plotParts( data, lbl, comp, domain=[0,1,0,1] ):
    ax = plt.subplot(111)
    plt.subplots_adjust(left=0.,bottom=0.25)
    t = sorted(data.keys())
    t_in = 0.
    dwis = list(set( [jj for (ii,jj) in data[min(t)].dw.keys() ] ))
    
    pl = plt.scatter([0],[0],s=10,linewidths=(0,0,0))
    plt.axis(domain)        
    
    axTime = plt.axes([0.25, 0.1, 0.65, 0.03]) 
    sTime = Slider(axTime, 'Time', min(t), max(t), valinit=min(t))
    
    
    def updateTime(val):
        t0 = sTime.val
        dt = abs(np.array(t)-t0)
        t_in = t[dt.argmin()]
        px = partVar(data,t_in,dwis,'Pos')
        pc = partVar(data,t_in,dwis,lbl)
        pl.set_offsets(px)  
        pl.set_array(pc)
        plt.draw()
            
    updateTime(t_in)
    sTime.on_changed(updateTime)
    
    plt.show()
    
    
#======================================================================
def partVar( data, t, dwis, lbl='Vel', comp='norm' ):
    if lbl=='Pos': 
        func = pPosition
    elif lbl=='Vel': 
        func = pVelocity
    elif lbl=='Stress': 
        func = pVonMises
    else:
        func = getPVar

    out = func( data, t, dwis[0], lbl, comp )
    for ii in range(1,len(dwis)):
        out = np.append(out, func(data,t,dwis[ii], lbl, comp), axis=0 )
        
    return out

#======================================================================
def getPVar( data, t, dwi, lbl, comp ):
    pv = copy( data[t].get(lbl,dwi) )
    if comp=='x': return pv[:,0]
    if comp=='y': return pv[:,1]
    if comp=='norm': return vnorm(pv)

#======================================================================
def pPosition( data, t, dwi, lbl=1, comp=1 ):
    px = copy( data[t].get('px',dwi) )
    return px

#======================================================================
def pVelocity( data, t, dwi, lbl=1, comp=1 ):
    pv = copy( data[t].get('pxI',dwi) )
    return vnorm(pv)

#======================================================================
def pVonMises( data, t, dwi, lbl=1, comp=1 ):
    pVS,pVol = data[t].getMult( ['pVS','pVol'], dwi )
    pS = [pVS[ii]/pVol[ii] for ii in range(len(pVol))]
    ms = np.array([vonMises(pp) for pp in pS])
    return ms

#======================================================================
def plotNodes( data, lbl, comp, dwi, domain=[0,1,0,1] ):
    ax = plt.subplot(111)
    plt.subplots_adjust(left=0.1,bottom=0.25)
    t = sorted(data.keys())
    
    pl = plt.scatter([0],[0],s=20,linewidths=(0,0,0),vmin=0,vmax=100)
    plt.axis(domain)        
    
    axTime = plt.axes([0.25, 0.1, 0.65, 0.03]) 
    sTime = Slider(axTime, 'Time', min(t), max(t), valinit=min(t))       
    
    def updateTime(val):
        t0 = sTime.val
        dt = abs(np.array(t)-t0)
        t_in = t[dt.argmin()]
        nx = nodeVar(data,t_in,dwi,'Pos')
        nc = nodeVar(data,t_in,dwi,lbl,comp)
        pl.set_offsets(nx)  
        pl.set_array(nc)
        plt.draw()
            
    updateTime(t[0])
    sTime.on_changed(updateTime)
    plt.show()
    
    
#======================================================================
def nodeVar( data, t, dwi, lbl='Vel', comp='norm' ):
    if lbl=='Pos': 
        func = nPosition
    elif lbl=='Vel': 
        func = nVelocity
    else:
        func = getNVar

    out = func( data, t, dwi, lbl, comp )        
    return out

#======================================================================
def getNVar( data, t, dwi, lbl, comp ):
    nv = copy( data[t].get(lbl,dwi) )
    if comp=='x': nvar = nv[:,0]
    if comp=='y': nvar = nv[:,1]
    if comp=='norm': nvar = vnorm(nv)
    return nvar

#======================================================================
def nPosition( data, t, dwi, lbl=1, comp=1 ):
    px = copy( data[t].get('gx',dwi) )
    return px

#======================================================================
def nVelocity( data, t, dwi, lbl=1, comp=1 ):
    pv = copy( data[t].get('gv',dwi) )
    return vnorm(pv)

#======================================================================    
def vnorm( x ):
    xx = np.sqrt((x*x).sum(axis=1)[:,np.newaxis])
    return np.reshape( xx, xx.size )

#======================================================================    
def vonMises( S ):
    return np.sqrt( S[0,0]*S[0,0] - S[0,0]*S[1,1] + S[1,1]*S[1,1] +
                    3.*S[1,0]*S[0,1] )