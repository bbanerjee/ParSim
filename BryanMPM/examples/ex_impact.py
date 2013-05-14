import numpy as np
import time
from mpm_imports import *
import sys


#===============================================================================
def bcZero( x ):
    return( np.zeros(2) )


#===============================================================================
def init(useCython):
    #  Initialize Simulation
    # Result File Name
    outputName = 'impact'
    outputDir = 'test_data/impact'

    
    # Domain Constants
    x0 = np.array([0.0,0.0])                 # Bottom left Corner
    x1 = np.array([3.,1.])                   # Top right corner
    nN = np.array([96,32])                   # Number of cells in each direction
    dx = (x1-x0)/nN                          # Cell size
    nG = 2                                   # Number of ghost nodes
    thick = 0.1                              # Domain thickness
    ppe = 2                                  # Particles per element 
    matid1 = 1;    matid2 = 2                # Material IDs
    initVel1 = np.array([15.,0.])            # Initial velocity vector #1
    initVel2 = np.array([0.,0.])             # Initial velocity vector #2
    gravity = np.array([0.,-9.8])            # Gravity Vector     
    

    # Material Properties
    E1 = 1.0e7;    nu1 = 0.3;    rho1 = 1.0e3;    
    vWave1 = np.sqrt( E1/rho1 )    
    matProps1 = { 'modulus':E1, 'poisson':nu1, 'density':rho1 }    
    matModelName = 'planeStrainNeoHookean'
        
        
    # Time Constants
    t0 = 0.0                                                  # Initial Time
    CFL = 0.4                                                 # CFL Condition
    dt = min(dx) * CFL / vWave1;                              # Time interval
    tf = 10.                                                  # Final time
    outputInterval = tf/400.                                  # Output interval


    # Create Data Warehouse, Patch, Shape functions
    dw = Dw( outputInterval, outputDir )
    patch = Patch( x0, x1, nN, nG, t0, tf, dt, thick, dw )
    shape = Shape( useCython )


    # Create boundary conditions
    bcVy = Bc( 'Y', 0.0, 'gv', bcZero )       # Zero velocity bc at y=0
    bcAy = Bc( 'Y', 0.0, 'ga', bcZero )       # Zero acceleration bc at y=0
    bcVx = Bc( 'X', 0.0, 'gv', bcZero )       # Zero velocity bc at y=0
    bcAx = Bc( 'X', 0.0, 'ga', bcZero )       # Zero acceleration bc at y=0
    patch.bcs = [bcVy,bcAy,bcVx,bcAx]           # BC list
 
 
    # Create Objects
    matList = []
    matList.append( Material( matid1, matProps1, matModelName, useCython ) )
    matList.append( Material( matid2, matProps1, matModelName, useCython ) )    
    
    dx0 = dx[0]
    density = matProps1['density']
    pt11 = np.array([1.2+2.0*dx0,0.3])        # Impacter bottom left corner
    pt12 = np.array([2.0+2.0*dx0,0.5])
    pt21 = np.array([2.0+4.0*dx0,0.0])        # Target bottom left corner
    pt22 = np.array([2.2+4.0*dx0,0.8])    
    geomutils.fillRectangle( pt11, pt12, ppe, patch, dw, matid1, density ) 
    geomutils.fillRectangle( pt21, pt22, ppe, patch, dw, matid2, density ) 
    
    
    # Initialize Data Warehouse and Object Velocities
    dw.initArrays(shape.nSupport)    
    mpm.updateMats( dw, patch, matList, shape )
    matList[0].setVelocity( dw, initVel1 )
    matList[1].setVelocity( dw, initVel2)    
    matList[0].setExternalAcceleration( dw, gravity )
    matList[1].setExternalAcceleration( dw, gravity )        
    
    print 'dt = ' + str(patch.dt)    
    return (dw, patch, matList, shape, outputName )


#===============================================================================
def stepTime( dw, patch, matList, shape, outputName ):
    # Advance through time
    
    tbegin = time.time()
    try:
        while( (patch.t < patch.tf) and patch.allInPatch(dw.px) ):
            dw.resetNodes()     
            mpm.timeAdvance( dw, patch, matList, shape )
            dw.saveDataAndAdvance( patch.dt, outputName )
    except JacobianError:
        print 'Negative Jacobian'
    
    tend = time.time()
    print (str(dw.idx) + ' iterations in: ' + readTime(tend-tbegin) 
           + ' t=' + str(patch.t) )

#===============================================================================            
def run(useCython=True):
    dw, patch, matList, shape, outputName = init(useCython)
    stepTime( dw, patch, matList, shape, outputName )
    return dw