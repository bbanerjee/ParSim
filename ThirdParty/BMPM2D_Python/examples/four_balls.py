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
# Initial Velocity Function
def initVel(x):
    dv = 1.
    if (x[0] > 0.12):
        dv = 0.
    return dv * 0.2 * np.array([1.,0.])


#===============================================================================
# Run the MPM simulation
#===============================================================================
def run( outputFile, saveDWs=False, useCython=True ):
    mpm = init(outputFile, useCython)
    mpmData = stepTime( mpm, saveDWs )
    return mpmData


#===============================================================================
# Initialize the simulation
#===============================================================================
def init( outputFile, useCython ):
    # Initialize Simulation
    # Result File Name + Directory
    outputDir = 'test_data/four_balls'
    
    #========================================
    # Domain Constants
    x0 = np.array([0.0,0.0]);                    # Bottom left corner
    x1 = np.array([0.25,0.25])                       # Top right corner
    nN = np.array([25,25])                       # Number of cells
    nG = 2                                       # Number of ghost nodes
    thick = 0.02                                  # Domain thickness
    ppe = 4                                      # Particles per element

    #========================================
    # Material Properties
    E1 = 1.24e6;    nu1 = 0.3;    rho1 = 8.0e3;    
    vWave1 = np.sqrt( E1/rho1 )
    props = {'modulus':E1, 'poisson':nu1, 'density':rho1 }
    matModel = 'planeStrainNeoHookean'

    #========================================    
    # Time Constants
    t0 = 0.0                                                # Initial Time
    CFL = 0.1                                               # CFL Condition
    dt = min((x1-x0)/nN) * CFL / vWave1;                    # Time interval
    tf = 2.                                                 # Final time
    outputInterval = 0.005                                  # Output interval  

    #========================================     
    # Create Data Warehouse, Patchs, Shape functions
    dw = Dw( outputDir, outputFile, outputInterval )
    patch = Patch( x0, x1, nN, nG, t0, tf, dt, thick, ppe )
    shape = Shape( useCython )

    #========================================
    # Create Objects
    mats = []    
    dwis = [1,2,3,4]
    r = 0.047134567;     
    xx0 = 0.06;  xx1 = 0.19;
    yy0 = 0.06;  yy1 = 0.19;
    x = [xx0,xx1+np.pi/1000,xx0,xx1+0.0025]
    y = [yy0,yy0,yy1,yy1]
    
    for ii in range(len(dwis)):
        cntr = np.array([x[ii],y[ii]])
        circ = gu.ellipseLvl( r, cntr )             
        dw.createGrid(dwis[ii],patch)
        mats.append(Material( props, matModel, dwis[ii], shape, useCython ))
        px, vol = gu.fillLvl( circ, patch )
        dw.addParticles( dwis[ii], px, vol, props['density'], shape.nSupport )

   
    #========================================
    # Create Contact
    mu = 0.5                                 # Friction Coefficient (if used)
    contacts = []
    contacts.append( FrictionContact([dwis[0],dwis[1]], patch, mu) )
    contacts.append( FrictionContact([dwis[2],dwis[3]], patch, mu) )


    #========================================
    # Create boundary conditions
    patch.bcs = []                                 #  BC list
    patch.bcs.append( Bc('Y', x0[1], 'gv', dwis, zero) )
    patch.bcs.append( Bc('Y', x1[1], 'gv', dwis, zero) )
    patch.bcs.append( Bc('X', x0[0], 'gv', dwis, zero) )
    patch.bcs.append( Bc('X', x1[0], 'gv', dwis, zero) )
    

    #========================================
    # Initialize Data Warehouse and Object Accelerations
    mpm2d.updateMats( dw, patch, mats )
    for mat in mats:
        mat.setVelocity( dw,  initVel )

    #========================================
    # Create dict 
    mpm = { 'dw': dw, 'patch': patch, 'mats': mats, 'contacts': contacts }

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

    #===========================================================================    
    # Check for Negative Jacobian Error
    ii = 0
    try:
        mpmData[dw.t] = copy(dw)          # Copy first dw
 
        #=======================================================================        
        # Loop while t<dt and all particles remain in the domain
        while( (patch.t < patch.tf) and inPatch ):
            dw.dumpData( patch.dt, mats )            
            mpm2d.timeAdvance( dw, patch, mats, contacts )
            
    except JacobianError:
        print 'Negative Jacobian'

    mpmData[dw.t] = copy(dw)            # Save last dw
 
    #===========================================================================
    # Print total time and number of iterations
    tend = time.time()    
    print (str(dw.idx) + ' iterations in: ' + readTime(tend-tbegin) 
            + ' t=' + str(patch.t) )
    return (mpmData)