import numpy as np
import time
from copy import deepcopy as copy
from itertools import izip, count
from mpm_imports import *
Shape = GIMP
Contact = FrictionContact

def zero(x): return 0.

#===============================================================================
# Run the MPM simulation
#===============================================================================
def run( outputFile, saveDWs=False, useCython=True ):
    mpm = init(outputFile, useCython)
    mpmData, t, pos = stepTime( mpm, saveDWs )
    return (mpmData,t,pos)


#===============================================================================
# Initialize the simulation
#===============================================================================
def init( outputFile, useCython ):
    # Initialize Simulation
    # Result File Name + Directory
    outputDir = 'test_data/slide'
    
    #========================================
    # Domain Constants
    x0 = np.array([0.0,0.0]);                    # Bottom left corner
    x1 = np.array([1.,1.])                       # Top right corner
    nN = np.array([20,20])                       # Number of cells
    nG = 2                                       # Number of ghost nodes
    thick = 0.1                                  # Domain thickness
    ppe = 4                                      # Particles per element

    #========================================
    # Material Properties
    E1 = 100.0e3;    nu1 = 0.3;    rho1 = 2.0e3;    
    vWave1 = np.sqrt( E1/rho1 )
    matProps1 = {'modulus':E1, 'poisson':nu1, 'density':rho1 }
    matProps = [matProps1]
    matModelName = 'planeStrainNeoHookean'

    #========================================    
    # Time Constants
    t0 = 0.0                                                  # Initial Time
    CFL = 0.2                                                 # CFL Condition
    dt = min((x1-x0)/nN) * CFL / vWave1;                      # Time interval
    tf = 1.                                                 # Final time
    outputInterval = dt                                    # Output interval  

    #========================================     
    # Create Data Warehouse, Patchs, Shape functions
    dw = Dw( outputDir, outputFile, outputInterval )
    patch = Patch( x0, x1, nN, nG, t0, tf, dt, thick, ppe )
    shape = Shape( useCython )

    #========================================
    # Create Objects
    mats = []    
    dwis = [1]
    r = 0.2
    cntr = np.array([0.2, 0.4])     # Cylinder center position
    circ,circNml = gu.ellipseLvl( r, cntr ) # Level set defining the cylinder and normal
    lvls = [circ]
    nmls = [circNml]
    for ii in range(len(dwis)):
        dens = matProps[ii]['density']
        dw.createGrid(dwis[ii],patch)
        mats.append(Material( matProps[ii], matModelName, dwis[ii], 
                              shape, useCython ))
        px, vol, nml = gu.fillLvl( lvls[ii], nmls[ii], patch )
        dw.addParticles( dwis[ii], px, vol, nml, dens, shape.nSupport )
   
    contacts = []

    #========================================
    # Create boundary conditions
    patch.bcs = []                                 #  BC list
    patch.bcs.append( Bc('Y', 0., 'gv', dwis, zero) )

    #========================================
    # Initialize Data Warehouse and Object Accelerations
    mpm2d.updateMats( dw, patch, mats )
    mats[0].setVelocity( dw,  initVel )

    #========================================
    # Find index of particle at center of cylinder
    cntrIdx = getCntr( dw, 1, cntr )

    #========================================
    # Create dict 
    mpm = { 'dw': dw, 'patch': patch, 'mats': mats, 'contacts': contacts,
            'cntrIdx': cntrIdx, 'cntrDwi': 1 }

    #========================================
    # Output time interval and return
    print 'dt = ' + str(patch.dt)        
    return mpm


#===============================================================================
# Step through time
#===============================================================================
def stepTime( mpm, saveDWs ):
    # Advance through time
    # if saveDWs, save the dw at each output step - otherwise just 1st + last
    dw = mpm['dw'];  patch = mpm['patch']; mats = mpm['mats']   # Load variables
    contacts = mpm['contacts']
    tbegin = time.time()
    mpmData = dict()
    pos = []; t = []    
    inPatch = True
    px0 = copy( findCntr( dw, mpm['cntrDwi'], mpm['cntrIdx'] ) )

    #========================================    
    # Check for Negative Jacobian Error
    try:
        ii=0
        if not saveDWs: mpmData[dw.t] = copy(dw)          # Copy first dw
 
        #========================================        
        # Loop while t<dt and all particles remain in the domain
        while( (patch.t < patch.tf) and inPatch ):
            mpm2d.timeAdvance( dw, patch, mats, contacts )
            
            if saveDWs and dw.checkSave(patch.dt): mpmData[daw.t] = copy(dw)
            dw.dumpData( patch.dt, mats )

            # Track cylinder center position
            t.append(dw.t)
            pos.append( np.linalg.norm(findCntr(dw,mpm['cntrDwi'],
                                                    mpm['cntrIdx'])-px0) )
            
            #========================================
            # Check that all particles are in the domain
            for mat in mats:
                if not patch.allInPatch(dw.get('px',mat.dwi)): inPatch = True
    except JacobianError:
        print 'Negative Jacobian'

    if not saveDWs: mpmData[dw.t] = copy(dw)            # Save last dw
 
    #========================================
    # Print total time and number of iterations
    tend = time.time()    
    print (str(dw.idx) + ' iterations in: ' + readTime(tend-tbegin) 
            + ' t=' + str(patch.t) )
    t = np.array(t)
    pos = np.array(pos)
    return (mpmData,t,pos)
    

#===============================================================================            
# Utility functions	
#===============================================================================
def plotParts( data, lbl, comp, domain=[0,1,0,1] ):
    plotutil.plotParts( data, lbl, comp, domain )
    
def plotNodes( data, lbl, comp, dwi, domain=[0,1,0,1] ):
    plotutil.plotNodes( data, lbl, comp, dwi, domain )    
    
def getCntr( dw, dwi, cntr ):
    pX = dw.get('pX', dwi)
    dist = vnorm(pX-cntr)
    return dist.argmin()

def findCntr( dw, dwi, cidx ):
    px = dw.get('px',dwi)
    return px[cidx]

def vnorm( x ):
    return np.sqrt((x*x).sum(axis=1)[:,np.newaxis])

def initVel(x):
    c = np.array([0.2,0.4])
    d = x-c
    r = np.linalg.norm(d)
    r = max(r,1.e-10)
    return 20.*np.array([d[1], -d[0]])+np.array([1.,0])


#===============================================================================
if __name__ =="__main__":
    run('test_')