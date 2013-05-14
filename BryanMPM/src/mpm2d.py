import time
import numpy as np

def timeAdvance( dw, patch, mats, sh ):
    # Advance timestep
    updateMats( dw, patch, mats, sh )
    applyExternalLoads( dw, patch, mats )
    interpolateParticlesToGrid( dw, patch, mats )
    computeStressTensor( dw, patch, mats )
    computeInternalForce( dw, patch, mats )
    computeAndIntegrateAcceleration( dw, patch, patch.tol )
    setGridBoundaryConditions( dw, patch, patch.tol )
    interpolateToParticlesAndUpdate( dw, patch, mats, sh )

    
def updateMats( dw, patch, mats, sh ):
    dw.resetNodes()
    for mat in mats:
        mat.getParticles(dw, patch, sh)
    
    
def applyExternalLoads( dw, patch, mats ):
    # Apply external loads to each material
    for mat in mats:
        mat.applyExternalLoads( dw, patch )

    
def interpolateParticlesToGrid( dw, patch, mats ):
    # Interpolate particle mass and momentum to the grid
    for mat in mats:
        mat.interpolateParticlesToGrid( dw, patch )    

    
def computeInternalForce( dw, patch, mats ):
    # Compute internal body forces
    for mat in mats:
        mat.computeInternalForce( dw, patch )

    
def computeAndIntegrateAcceleration( dw, patch, tol ):
    # Integrate grid acceleration
    dt = patch.dt
    a_leap = 1                                        # Initializes leap-frog
    if( patch.it == 0 ): a_leap = 0.5

    gm = dw.getData( 'gm' )                           # Mass
    gw = dw.getData( 'gw' )                           # Momentum
    gfi = dw.getData( 'gfi' )                         # Internal Force
    gfe = dw.getData( 'gfe' )                         # External Force
    gv = dw.getData( 'gv' )                           # Velocity
    ga = dw.getData( 'ga' )
    
    gm[:] += tol
    gv[:] = gw/gm
    ga[:] = a_leap * (gfe+gfi)/gm
    gv[:] += ga*dt
    

def setGridBoundaryConditions( dw, patch, tol ):
    for bc in patch.bcs:
        bc.setBoundCond( dw, patch, tol )


def computeStressTensor( dw, patch, mats ):
    for mat in mats:
        mat.computeStressTensor( dw, patch )

    
def interpolateToParticlesAndUpdate( dw, patch, mats, sh ):
    for mat in mats:
        mat.interpolateToParticlesAndUpdate( dw, patch )
        
    patch.stepTime()