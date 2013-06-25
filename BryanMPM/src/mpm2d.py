import time
import numpy as np

def timeAdvance( dw, patch, mats, contacts=[] ):
    # Advance timestep
    updateMats( dw, patch, mats )
    applyExternalLoads( dw, patch, mats )
    interpolateParticlesToGrid( dw, patch, mats )
    exchMomentumInterpolated( dw, contacts )
    computeStressTensor( dw, patch, mats )
    computeInternalForce( dw, patch, mats )
    exchForceInterpolated( dw, contacts )    
    computeAndIntegrateAcceleration( dw, patch, mats )
    exchMomentumIntegrated( dw, contacts )
    setGridBoundaryConditions( dw, patch )
    interpolateToParticlesAndUpdate( dw, patch, mats )

    
def updateMats( dw, patch, mats ):
    for mat in mats:
        mat.updateContributions(dw, patch)
    
    
def applyExternalLoads( dw, patch, mats ):
    # Apply external loads to each material
    for mat in mats:
        mat.applyExternalLoads( dw, patch )

    
def interpolateParticlesToGrid( dw, patch, mats ):
    # Interpolate particle mass and momentum to the grid
    for mat in mats:
        mat.interpolateParticlesToGrid( dw, patch )    


def exchMomentumInterpolated( dw, contacts ):
    # Exchange Interpolated Momentum
    for contact in contacts:
        contact.exchMomentumInterpolated( dw )

def computeStressTensor( dw, patch, mats ):
    # Compute Material Stress Tensors
    for mat in mats:
        mat.computeStressTensor( dw, patch )
    
def computeInternalForce( dw, patch, mats ):
    # Compute internal body forces
    for mat in mats:
        mat.computeInternalForce( dw, patch )

def exchForceInterpolated( dw, contacts ):
    # Exchange Interpolated Momentum
    for contact in contacts:
        contact.exchForceInterpolated( dw )
    
def computeAndIntegrateAcceleration( dw, patch, mats ):
    # Integrate grid acceleration
    for mat in mats:
        mat.computeAndIntegrateAcceleration( dw, patch, patch.tol )    

def exchMomentumIntegrated( dw, contacts ):
    # Exchange Interpolated Momentum
    for contact in contacts:
        contact.exchMomentumInterpolated( dw )

def setGridBoundaryConditions( dw, patch ):
    # Set boundary conditions
    for bc in patch.bcs:
        bc.setBoundCond( dw, patch, patch.tol )
    
def interpolateToParticlesAndUpdate( dw, patch, mats ):
    # Interpolate velocity/accel/deformation vals to Particles and Move
    for mat in mats:
        mat.interpolateToParticlesAndUpdate( dw, patch )
        
    patch.stepTime()