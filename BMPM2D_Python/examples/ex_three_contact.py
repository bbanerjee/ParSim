import numpy as np
import time
from copy import deepcopy as copy
from itertools import izip, count
from mpm_imports import *
Shape = Quad
Contact = FrictionContact

#===============================================================================
# Initial Velocity Function
def initVel(x):
    dv = 1.
    if (x[1]+x[0] > 2.):
        dv = -1.
    return dv * 0.1 * np.array([1.,1.])

def vel0(x): return  0.1 * np.array([1.,1.])
def vel1(x): return  0.0 * np.array([1.,1.])
def vel2(x): return -0.1 * np.array([1.,1.])


#===============================================================================
# Initialize the simulation
def init( outputName, useCython ):
    # Initialize Simulation
    # Result File Name + Directory
    outputDir = 'test_data/three_contact'

    
    # Domain Constants
    x0 = np.array([0.0,0.0]);                    # Bottom left corner
    x1 = np.array([5.,5.])                       # Top right corner
    nN = np.array([100,100])                       # Number of cells
    nG = 2                                       # Number of ghost nodes
    thick = 0.1                                  # Domain thickness
    ppe = 4                                      # Particles per element
    mu = 0.3                                     # Friction Coefficient

    # Material Properties
    E1 = 100.0e3;    nu1 = 0.3;    rho1 = 1.0e3;    
    vWave1 = np.sqrt( E1/rho1 )
    matProps1 = {'modulus':E1, 'poisson':nu1, 'density':rho1 }
    matProps2 = {'modulus':E1*10, 'poisson':nu1, 'density':rho1*10 }
    matProps = [matProps1, matProps1, matProps2]
    matModelName = 'planeStrainNeoHookean'

    
    # Time Constants
    t0 = 0.0                                                  # Initial Time
    CFL = 0.4                                                 # CFL Condition
    dt = min((x1-x0)/nN) * CFL / vWave1;                      # Time interval
    tf = 10.                                                  # Final time
    outputInterval = 0.01                                     # Output interval
    
     
    # Create Data Warehouse, Patchs, Shape functions
    dw = Dw( outputDir, outputName, outputInterval )
    patch = Patch( x0, x1, nN, nG, t0, tf, dt, thick, ppe )
    shape = Shape( useCython )


    # Create boundary conditions
    patch.bcs = []                                 #  BC list
    
    # Create Objects
    mats = []    
    dwis = [1,2,3]
    ctrs = [0.75, 1.5, 2.25]
    
    for ii in range(len(dwis)):
        dens = matProps[ii]['density']
        dw.createGrid(dwis[ii],patch)
        mats.append(Material( matProps[ii], matModelName, dwis[ii], 
                              shape, useCython ))
        circ = gu.ellipseLvl( 0.2, np.array([ctrs[ii],ctrs[ii]]) )       # r, center
        px, vol = gu.fillLvl( circ, patch )
        dw.addParticles( dwis[ii], px, vol, dens, shape.nSupport )
    
    # Create Contact
    contacts = []
    contacts.append( Contact([dwis[0],dwis[1]], mu, dt, ppe) )
    contacts.append( Contact([dwis[1],dwis[2]], mu, dt, ppe) )

    # Initialize Data Warehouse and Object Velocities
    mpm.updateMats( dw, patch, mats )
    vels = [1., 0., 0.]    
    for (mat,vel) in izip(mats,vels):
        mat.setVelocity( dw,  vel*np.array([1.,1.]) )
    
    print 'dt = ' + str(patch.dt)        
    return (dw, patch, mats, contacts )


#===============================================================================
def stepTime( dw, patch, mats, contacts, useSave ):
    # Advance through time
    tbegin = time.time()
    mpmData = dict()
    try:
        if not useSave: mpmData[dw.t] = copy(dw)
        inPatch = True
        while( (patch.t < patch.tf) and inPatch ):
            mpm.timeAdvance( dw, patch, mats, contacts )
            if useSave and dw.checkSave(patch.dt): mpmData[dw.t] = copy(dw)
            
            dw.saveData( patch.dt, mats )

            for mat in mats:
                if not patch.allInPatch(dw.get('px',mat.dwi)): inPatch = False
    except JacobianError:
        print 'Negative Jacobian'

    if not useSave: mpmData[dw.t] = copy(dw)
    
    tend = time.time()
    print (str(dw.idx) + ' iterations in: ' + readTime(tend-tbegin) 
            + ' t=' + str(patch.t) )
            
    return mpmData
    

#===============================================================================            
def run( output='two', useSave=True, useCython=True ):
    dw, patch, mats, contacts = init( output, useCython )
    mpmData = stepTime( dw, patch, mats, contacts, useSave )
    return mpmData

def plot( data, domain=[0,1,0,1] ):
    plotutil.plot( data, domain )