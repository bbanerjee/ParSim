#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
#

import numpy as np
import numpy.linalg as la
import time
from copy import deepcopy as copy
from itertools import izip, count
from mpm_imports import *
Shape = GIMP

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
    outputDir = 'test_data/bounce'
    
    #========================================
    # Domain Constants
    x0 = np.array([0.0,0.0]);                    # Bottom left corner
    x1 = np.array([0.65,0.65])                       # Top right corner
    nN = np.array([65,65])                       # Number of cells
    nG = 2                                       # Number of ghost nodes
    thick = 0.02                                  # Domain thickness
    ppe = 4                                      # Particles per element

    #========================================
    # Material Properties
    E1 = 1.24e6;    nu1 = 0.3;    rho1 = 8.0e3;    
    vWave1 = np.sqrt( E1/rho1 )
    matProps1 = {'modulus':E1, 'poisson':nu1, 'density':rho1 }
    matProps2 = {'modulus':E1*1000, 'poisson':nu1, 'density':rho1*100000 }
    matProps = [matProps1, matProps1, matProps1, matProps1, matProps1, matProps2]
    matModelName = 'planeStrainNeoHookean'

    #========================================    
    # Time Constants
    t0 = 0.0                                                  # Initial Time
    CFL = 0.15                                                 # CFL Condition
    dt = min((x1-x0)/nN) * CFL / vWave1;                      # Time interval
    tf = 2.                                                 # Final time
    outputInterval = 0.002                                  # Output interval  

    #========================================     
    # Create Data Warehouse, Patchs, Shape functions
    dw = Dw( outputDir, outputFile, outputInterval )
    patch = Patch( x0, x1, nN, nG, t0, tf, dt, thick, ppe )
    shape = Shape( useCython )

    #========================================
    # Create Objects
    mats = []    
    dwis = [1,2,3,4,5,6]
    y0 = 0.02;  r = 0.04;  offset = 0.1     # Plane height, radius, vert offset
    cntr1 = np.array([2.*r, y0+r+offset])      # Cylinder center position
    cntr2 = np.array([5.*r, y0+r+offset])      # Cylinder center position
    cntr3 = np.array([8.*r, y0+r+offset])      # Cylinder center position
    cntr4 = np.array([11.*r, y0+r+offset])      # Cylinder center position    
    cntr5 = np.array([14.*r, y0+r+offset])      # Cylinder center position        
    circ1 = gu.ellipseLvl( r, cntr1 )           # Level set defining the cylinder and normal
    circ2 = gu.ellipseLvl( r, cntr2 )           # Level set defining the cylinder and normal
    circ3 = gu.ellipseLvl( r, cntr3 )           # Level set defining the cylinder and normal
    circ4 = gu.ellipseLvl( r, cntr4 )           # Level set defining the cylinder and normal
    circ5 = gu.ellipseLvl( r, cntr5 )           # Level set defining the cylinder and normal
    def plane(x):  return y0 - x[1]            # Level set defining the plane    
    def planeNml(x): return np.array([0.,1.])
    lvls = [circ1, circ2, circ3, circ4, circ5, plane]
    for ii in range(len(dwis)):
        dens = matProps[ii]['density']
        dw.createGrid(dwis[ii],patch)
        mats.append(Material( matProps[ii], matModelName, dwis[ii], 
                              shape, useCython ))
        px, vol = gu.fillLvl( lvls[ii], patch )
        dw.addParticles( dwis[ii], px, vol, dens, shape.nSupport )
   
    #========================================
    # Create Contact
    mu = 0.5                                 # Friction Coefficient (if used)
    contacts = []
    contacts.append( FrictionContact([dwis[5],dwis[0]], patch, mu) )
    contacts.append( FrictionlessContact([dwis[5],dwis[1]], patch) )
    contacts.append( FrictionContactTest([dwis[5],dwis[2]], patch, mu) )
    contacts.append( FreeContact([dwis[5],dwis[3]], patch) )
    contacts.append( FreeContact([dwis[5],dwis[4]], patch) )

    #========================================
    # Create boundary conditions
    patch.bcs = []                                 #  BC list
    patch.bcs.append( Bc('Y', 0., 'gv', dwis, zero) )
    patch.bcs.append( Bc('X', 4., 'gv', dwis, zero) )
    

    #========================================
    # Initialize Data Warehouse and Object Accelerations
    mpm2d.updateMats( dw, patch, mats )
    g_angle = 0. * np.pi/180.
    g = 9.8 * np.array( [np.sin(g_angle), -np.cos(g_angle)] )
    accs = [g, g, g, g, g, 0.]
    for (mat,acc) in izip(mats,accs):
        mat.setExternalAcceleration( dw,  acc )        
    
    #========================================
    # Find index of particle at center of cylinder
    cntrIdx = getCntr( dw, 1, cntr1 )

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
    ii = 0
    try:
        if not saveDWs: mpmData[dw.t] = copy(dw)          # Copy first dw
 
        #========================================        
        # Loop while t<dt and all particles remain in the domain
        while( (patch.t < patch.tf) and inPatch ):
            if saveDWs and dw.checkSave(patch.dt): mpmData[daw.t] = copy(dw)
            dw.dumpData( patch.dt, mats )            

            # Track cylinder center position
            t.append(dw.t)
            pos.append(la.norm(findCntr(dw,mpm['cntrDwi'], mpm['cntrIdx'])-px0))
            mpm2d.timeAdvance( dw, patch, mats, contacts )
            
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


#===============================================================================
if __name__ =="__main__":
    run('test_')
