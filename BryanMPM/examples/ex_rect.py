import numpy as np
import time
from mpm_imports import *
import copy

#===============================================================================
def init():
    #  Initialize Simulation
    # Result File Name
    fName = 'rect'
    fDir = 'test_data/rect'
    
    # Domain Constants
    x0 = np.array([0.0,0.0])                 #  Bottom left Corner
    x1 = np.array([1.0,1.0])                 #  Top right corner
    nN = np.array([10,10])                   #  Number of cells in each direction
    dx = (x1-x0)/nN
    nG = 2                                   #  Number of ghost nodes
    thick = 0.1                              #  Domain thickness
    ppe = 2                                  #  Particles per element 

    # Material Properties
    mProps1 = {'modulus':1.0e3, 'poisson':0.3, 'density':1.0e3 }
    
    modelName = 'planeStrainNeoHookean'
    mat1 = Material( matid=1, props=mProps1, model=modelName, exload=0)
    mats = [mat1]
    
    # Time Constants
    vw = np.sqrt( mProps1['modulus']/mProps1['density'] )   # Wave speed
    t0 = 0.0;    CFL = 0.4
    dt = min(dx) * CFL / vw
    tf = 1000*dt
    
    # Create Data Warehouse, Patch, Shape functions
    dw = Dw( t=0.0, idx=0, sidx=0, tout=t_out, ddir=fDir )
    pch = Patch( x0, x1, nN, nG, t0, tf, dt, thick, dw )
    sh = Shape()

    # Create boundary conditions
    pch.bcs = []                           #  BC list
    
    # Create Rectangle
    pt1 = np.array([0.4,0.4]);  pt2 = np.array([0.6,0.6])  # Corners
    matid1 = 1                                             # Material ID
    geomutils.fillRectangle( pt1, pt2, ppe, pch, dw, 1, mProps1['density'])
    
    mpm.updateMats( dw, pch, mats, sh )
    v0 = np.array([0.1,0.0])                               # Initial velocity
    mat1.setVelocity( dw, v0 )
    
    dw.initArrays()
    
    return (dw, pch, mats, sh, fName )


#===============================================================================
def stepTime( dw, patch, mats, shape, saveName ):
    # Advance through time
    dws = []
    while( (patch.t < patch.tf) and patch.allInPatch(dw.px) ):
        dw.resetNodes()
        t_out = patch.dt     
        mpm.timeAdvance( dw, patch, mats, shape )
        dw.saveDataAndAdvance( patch.dt, saveName )
        dws.append(copy.deepcopy(dw))
    return dws

#===============================================================================            
def run():
    dw, pch, mats, sh, fname = init()
    dws = stepTime( dw, pch, mats, sh, fname )
    return dws