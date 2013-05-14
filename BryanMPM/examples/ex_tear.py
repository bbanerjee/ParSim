import numpy as np
import time
from mpm_imports import *
import sys


#===============================================================================
def bcZero( x ):
    return( np.zeros(2) )


#===============================================================================
def init(bCython):
    if bCython: Shape = Shape_c
    
    #  Initialize Simulation
    # Result File Name
    fName = 'tear'
    fDir = 'test_data/tear'
    
    # Domain Constants
    x0 = np.array([0.0,0.0])                 #  Bottom left Corner
    x1 = np.array([3.,1.])                   #  Top right corner
    nN = np.array([96,32])                   #  Number of cells in each direction
    dx = (x1-x0)/nN
    nG = 2                                   #  Number of ghost nodes
    thick = 0.1                              #  Domain thickness
    ppe = 2                                  #  Particles per element 

    # Material Properties
    mProps1 = {'modulus':1.0e7, 'poisson':0.3, 'density':100.0e3, 
               'maxStress': 1.5e6}
    mProps2 = {'modulus':1.0e7, 'poisson':0.3, 'density':1.0e3, 
                   'maxStress': 1.5e6}    
    vw = np.sqrt( mProps1['modulus']/mProps2['density'] )   #  Wave speed
    modelName = 'planeStrainNeoHookeanMaxStress'
    mat1 = Material( matid=1, props=mProps1, model=modelName, exload=0,
                     bCython=bCython, bIgnoreNeg=True )
    mat2 = Material( matid=2, props=mProps2, model=modelName, exload=0,
                     bCython=bCython, bIgnoreNeg=True )
    mat3 = Material( matid=3, props=mProps2, model=modelName, exload=0,
                     bCython=bCython, bIgnoreNeg=True )
    mats = [mat1,mat2,mat3]
     
    # Time Constants
    t0 = 0.0
    CFL = 0.4
    dt = min(dx) * CFL / vw
    tf = 0.15
    t_out = tf/100.            
    
    # Create Data Warehouse
    dw = Dw( t=0.0, idx=0, sidx=0, save_interval = t_out, ddir=fDir )
    
    # Create Patch
    pch = Patch( x0, x1, nN, nG, t0, tf, dt, thick, dw )
    
    # Create shape function
    sh = Shape()

    # Create boundary conditions
    bcVy = Bc( 'Y', 0.0, 'gv', bcZero )       # Zero velocity bc at y=0
    bcAy = Bc( 'Y', 0.0, 'ga', bcZero )       # Zero acceleration bc at y=0
    bcVx = Bc( 'X', 0.0, 'gv', bcZero )       # Zero velocity bc at x=0
    bcAx = Bc( 'X', 0.0, 'ga', bcZero )       # Zero acceleration bc at x=0
    pch.bcs = [bcVy,bcAy,bcVx,bcAx]           # BC list
    pch.bcs = []
 
    # Create Rectangles
    dx0 = dx[0]
    pt11 = np.array([0.2+2.0*dx0,0.3])
    pt12 = np.array([0.8+3.0*dx0,0.5])
    pt21 = np.array([0.8+3.0*dx0,0.3])
    pt22 = np.array([1.0+2.0*dx0,0.5])    
    pt31 = np.array([1.0+4.0*dx0,0.1])
    pt32 = np.array([1.2+4.0*dx0,0.9])    
    matid1 = 1
    matid2 = 2
    matid3 = 3
    geomutils.fillRectangle( 
        pt11, pt12, ppe, pch, dw, matid1, mProps1['density'])     # Impacter
    geomutils.fillRectangle( 
        pt21, pt22, ppe, pch, dw, matid2, mProps2['density'])     # Impacter head
    geomutils.fillRectangle( 
        pt31, pt32, ppe, pch, dw, matid3, mProps2['density'])     # Target    
    dw.initArrays(sh.nSupport)
    
    v0 = np.array([15.0,0.0])
    v1 = np.array([0.0,0.0])
    gravity = np.array([0.0,-9.8])
    
    
    mpm.updateMats( dw, pch, mats, sh )
    mat1.setVelocity( dw, v0 )
    mat2.setVelocity( dw, v1 )
    mat3.setVelocity( dw, v1 )
    mat1.setExternalLoad( dw, gravity )
    mat2.setExternalLoad( dw, gravity )
    mat3.setExternalLoad( dw, gravity )
    
    
    print 'dt = ' + str(pch.dt)    
    return (dw, pch, mats, sh, fName )


#===============================================================================
def stepTime( dw, patch, mats, shape, saveName ):
    # Advance through time
    tbegin = time.time()
    try:
        while( (patch.t < patch.tf) and patch.allInPatch(dw.px) ):
            dw.resetNodes()     
            mpm.timeAdvance( dw, patch, mats, shape )
            dw.saveDataAndAdvance( patch.dt, saveName )
    except JacobianError:
        print 'Negative Jacobian'
    
    tend = time.time()
    print (str(dw.idx) + ' iterations in: ' + readTime(tend-tbegin) 
           + ' t=' + str(patch.t) )

#===============================================================================            
def run(bCython):
    dw, pch, mats, sh, fname = init(bCython)
    stepTime( dw, pch, mats, sh, fname )
    return dw