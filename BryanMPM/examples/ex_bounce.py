import numpy as np
import time
from mpm_imports import *

#===============================================================================
def bcZero( x ):
    return( np.zeros(2) )

#===============================================================================
def init():
    # Initialize Simulation
    # Result File Name + Directory
    fName = 'bounce'
    fDir = 'test_data/bounce_data'
    
    # Domain Constants
    x0 = np.array([0.0,0.0]);                # Bottom left corner
    x1 = np.array([1.,1.])                   # Top right corner
    nN = np.array([20,20])                   # Number of cells
    dx = (x1-x0)/nN                          # Cell size
    nG = 2                                   # Number of ghost nodes
    thick = 0.1                              # Domain thickness
    ppe = 2                                  # Particles per element 

    # Material Properties
    mProps1 = {'modulus':1.0e3, 'poisson':0.3, 'density':1.0e3 }
    modelName = 'planeStrainNeoHookean'
    mat1 = Material( matid=1, props=mProps1, model=modelName, exload=0)
    mat2 = Material( matid=2, props=mProps1, model=modelName, exload=0)
    mats = [mat1,mat2]
    
    # Time Constants
    vw = np.sqrt( mProps1['modulus']/mProps1['density'] )   # Wave speed    
    t0 = 0.0;    CFL = 0.4
    dt = min(dx) * CFL / vw;
    tf = 10;  dt_save = 0.1      
    
    # Create Data Warehouse, Patch, Shape functions
    dw = Dw( t=0.0, idx=0, sidx=0, tout=dt_save, ddir=fDir )
    pch = Patch( x0, x1, nN, nG, t0, tf, dt, thick, dw )
    sh = Shape()

    # Create boundary conditions
    bcx1 = Bc( 'X', 0.0, 'gv', bcZero )       # Zero velocity bc at x=0
    bcx2 = Bc( 'X', 1.0, 'gv', bcZero )       # Zero velocity bc at x=1
    bcy1 = Bc( 'Y', 0.0, 'gv', bcZero )       # Zero velocity bc at x=0
    bcy2 = Bc( 'Y', 1.0, 'gv', bcZero )       # Zero velocity bc at x=1
    pch.bcs = [bcx1,bcx2,bcy1,bcy2]           # BC list                        #  BC list
    
    # Create Circles
    pt1 = np.array([0.25,0.25])                  # Center 1
    pt2 = np.array([0.75,0.75])                  # Center 2
    r = np.array([0.0,0.2])                      # Radius (inner and outer)
    matid1 = 1;  matid2 = 2                      # Material IDs
    geomutils.fillAnnulus( pt1,r[0],r[1],ppe,pch,dw,matid1,mProps1['density'] )
    geomutils.fillAnnulus( pt2,r[0],r[1],ppe,pch,dw,matid2,mProps1['density'] )
    
    mpm.updateMats( dw, pch, mats, sh )
    v0 = np.array([0.1,0.1])                     # Initial Velocity    
    mat1.setVelocity( dw, v0 )
    mat2.setVelocity( dw, -v0 )

    dw.initArrays()
    
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
    except JacobError:
        print 'Negative Jacobian'
            
    tend = time.time()
    print (str(dw.idx) + ' iterations in: ' + readTime(tend-tbegin) 
            + ' t=' + str(patch.t) )
    

#===============================================================================            
def run():
    dw, pch, mats, sh, fname = init()
    stepTime( dw, pch, mats, sh, fname )
    return dw