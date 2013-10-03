import numpy as np
import time
from copy import deepcopy as copy
from mpm_imports import *
Shape = GIMP

#===============================================================================
# Initial Velocity Function
def initVel(x):
    dv = 1.
    if (x[1]+x[0] > 2.):
        dv = -0.
    return dv * 0.1 * np.array([1.,1.])


#===============================================================================
# Initialize the simulation
def init( outputName, useCython ):
    # Initialize Simulation
    # Result File Name + Directory
    outputDir = 'test_data/two_contact'

    
    # Domain Constants
    x0 = np.array([0.0,0.0]);                    # Bottom left corner
    x1 = np.array([2.,2.])                       # Top right corner
    nN = np.array([40,40])                       # Number of cells
    nG = 2                                       # Number of ghost nodes
    thick = 1.                                  # Domain thickness
    ppe = 4                                      # Particles per element


    # Material Properties
    E1 = 1.0e3;    nu1 = 0.3;    rho1 = 1.0e3;    
    vWave1 = np.sqrt( E1/rho1 )
    matProps1 = {'modulus':E1, 'poisson':nu1, 'density':rho1 }
    matModelName = 'planeStrainNeoHookean'

    
    # Time Constants
    t0 = 0.0                                                  # Initial Time
    CFL = 0.4                                                 # CFL Condition
    dt = min((x1-x0)/nN) * CFL / vWave1;                      # Time interval
    tf = 8.                                                  # Final time
    outputInterval = 0.05                                     # Output interval
    
     
    # Create Data Warehouse, Patchs, Shape functions
    dw = Dw( outputDir, outputName, outputInterval )
    patch = Patch( x0, x1, nN, nG, t0, tf, dt, thick, ppe )
    shape = Shape( useCython )


    # Create boundary conditions
    patch.bcs = []                                 #  BC list
    
    # Create Objects
    mats = []    
    dwis = [1]    
    dw.createGrid(dwis[0],patch)
    mats.append(Material( matProps1, matModelName, dwis[0], shape, useCython ))
    
    density = matProps1['density']
    circ1 = gu.ellipseLvl( 0.2, np.array([0.75,0.75]) )       # r, center
    circ2 = gu.ellipseLvl( 0.2, np.array([1.75,1.75]) )    
    #circ1 = gu.rectLvl( np.array([0.55,0.58]), np.array([0.95,0.98]) ) # x0, x1
    #circ2 = gu.rectLvl( np.array([1.05,1.05]), np.array([1.45,1.45]) )
    px1, vol1 = gu.fillLvl( circ1, patch )
    px2, vol2 = gu.fillLvl( circ2, patch )    
    dw.addParticles( dwis[0], px1, vol1, density, shape.nSupport )
    dw.addParticles( dwis[0], px2, vol2, density, shape.nSupport )
    
    # Create Contact
    contacts = []

    # Initialize Data Warehouse and Object Velocities
    mpm.updateMats( dw, patch, mats )
    mats[0].setVelocity( dw,  initVel )
  

    print 'dt = ' + str(patch.dt)        
    return (dw, patch, mats, contacts )


#===============================================================================
def stepTime( dw, patch, mats, contacts ):
    # Advance through time
    tbegin = time.time()
    mpmData = dict()
    try:
        inPatch = True
        while( (patch.t < patch.tf) and inPatch ):
            mpm.timeAdvance( dw, patch, mats, contacts )
            if dw.checkSave(patch.dt): mpmData[dw.t] = copy(dw)
            dw.saveData( patch.dt, mats )
            for mat in mats:
                if not patch.allInPatch(dw.get('px',mat.dwi)): inPatch = False
    except JacobianError:
        print 'Negative Jacobian'
            
    tend = time.time()
    print (str(dw.idx) + ' iterations in: ' + readTime(tend-tbegin) 
            + ' t=' + str(patch.t) )
            
    return mpmData
    

#===============================================================================            
def run( output='two', useCython=True ):
    dw, patch, mats, contacts = init( output, useCython )
    mpmData = stepTime( dw, patch, mats, contacts )
    return mpmData

def plot( data, domain=[0,1,0,1] ):
    plotutil.plot( data, domain )