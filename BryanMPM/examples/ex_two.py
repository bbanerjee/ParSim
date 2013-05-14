import numpy as np
import time
from mpm_imports import *

#===============================================================================
def init( useCython ):
    # Initialize Simulation
    # Result File Name + Directory
    outputName = 'two'
    outputDir = 'test_data/two'

    
    # Domain Constants
    x0 = np.array([0.0,0.0]);                    # Bottom left corner
    x1 = np.array([2.,2.])                       # Top right corner
    nN = np.array([40,40])                       # Number of cells
    dx = (x1-x0)/nN                              # Cell size
    nG = 2                                       # Number of ghost nodes
    thick = 0.1                                  # Domain thickness
    ppe = 4                                      # Particles per element
    initVelocity = 0.1 * np.array([1.,1.])       # Initial velocity vector
    matid1 = 1;    matid2 = 2                    # Material IDs


    # Material Properties
    E1 = 1.0e3;    nu1 = 0.3;    rho1 = 1.0e3;    
    vWave1 = np.sqrt( E1/rho1 )
    matProps1 = {'modulus':E1, 'poisson':nu1, 'density':rho1 }
    matModelName = 'planeStrainNeoHookean'

    
    # Time Constants
    t0 = 0.0                                                  # Initial Time
    CFL = 0.4                                                 # CFL Condition
    dt = min(dx) * CFL / vWave1;                              # Time interval
    tf = 10.                                                  # Final time
    outputInterval = 0.05                                     # Output interval
    
     
    # Create Data Warehouse, Patchs, Shape functions
    dw = Dw( outputInterval, outputDir )
    patch = Patch( x0, x1, nN, nG, t0, tf, dt, thick, dw )
    shape = Shape( useCython )


    # Create boundary conditions
    patch.bcs = []                                 #  BC list

    
    # Create Objects
    matList = []
    matList.append( Material( matid1, matProps1, matModelName, useCython ) )
    matList.append( Material( matid2, matProps1, matModelName, useCython ) )    
    
    center1 = np.array([0.75,0.75])
    center2 = np.array([1.25,1.25])
    radii = np.array([0.0,0.2])
    density = matProps1['density']
    geomutils.fillAnnulus( center1, radii, ppe, patch, dw, matid1, density )
    geomutils.fillAnnulus( center2, radii, ppe, patch, dw, matid2, density )


    # Initialize Data Warehouse and Object Velocities
    dw.initArrays(shape.nSupport)    
    mpm.updateMats( dw, patch, matList, shape )
    matList[0].setVelocity( dw,  initVelocity )
    matList[1].setVelocity( dw, -initVelocity )

    print 'dt = ' + str(patch.dt)        
    return (dw, patch, matList, shape, outputName )


#===============================================================================
def stepTime( dw, patch, mats, shape, outputName ):
    # Advance through time
    tbegin = time.time()
    try:
        while( (patch.t < patch.tf) and patch.allInPatch(dw.px) ):
            mpm.timeAdvance( dw, patch, mats, shape )
            dw.saveDataAndAdvance( patch.dt, outputName )
    except JacobianError:
        print 'Negative Jacobian'
            
    tend = time.time()
    print (str(dw.idx) + ' iterations in: ' + readTime(tend-tbegin) 
            + ' t=' + str(patch.t) )
    

#===============================================================================            
def run( useCython=True ):
    dw, patch, matList, shape, outputName = init( useCython )
    stepTime( dw, patch, matList, shape, outputName )
    return dw