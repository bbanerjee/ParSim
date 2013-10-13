/*
 * MPMMaterial.cc
 *
 *  Created on: 11/10/2013
 *      Author: banerjee
 */

#include "MPMMaterial.h"

using namespace BrMPM;

// Material - holds update functions - default is deformable
// overridden by RigidMaterial for rigid materials
MPMMaterial::MPMMaterial() {
    self.props = props
    self.dwi = dwi
    self.shape = shape

    try:
        self.ignoreNegJ = props['ignoreNegJ']
    except Exception:
        self.ignoreNegJ = False

    if useCython:
        self.util = util_c
        self.mmodel = mmodel_c
    else:
        self.util = util
        self.mmodel = mmodel

    self.mm = self.mmodel.MaterialModel( model, props )

}

MPMMaterial::~MPMMaterial() {
  // TODO Auto-generated destructor stub
}


    void BrMPM::MPMMaterial::updateContributions(MPMDatawarehouseP& dw, MPMPatchP& patch)
    {
    dw.zeroGrid( self.dwi )
    self.shape.updateContribList( dw, patch, self.dwi )
    }

    void BrMPM::MPMMaterial::setVelocity(MPMDatawarehouseP& dw,
        MPMParticleVar<Vector>& v) {
    pw,pm,px = dw.getMult( ['pw','pm','px'], self.dwi )

    for (ii,pxi,pmi) in izip(count(),px,pm):
        if isfunction(v):
            pw[ii] = v(pxi) * pmi
        else:
            pw[ii] = v * pmi
    }

    void BrMPM::MPMMaterial::setExternalLoad(MPMDatawarehouseP& dw,
        MPMParticleVar<Vector> fe) {
    pfe = dw.get( 'pfe', self.dwi )
    for pfei in pfe:
        pfei = fe
    }

    void BrMPM::MPMMaterial::setExternalAcceleration(MPMDatawarehouseP& dw,
        MPMParticleVar<Vector> acc) {
    pfe,pm = dw.getMult( ['pfe','pm'], self.dwi )
    for (ii,pmi) in izip(count(),pm):
        pfe[ii] = acc*pmi
    }

    void BrMPM::MPMMaterial::applyExternalLoads(MPMDatawarehouseP& dw,
        MPMPatchP& patch) {
    # Apply external loads to each material
    cIdx,cW = dw.getMult( ['cIdx','cW'], self.dwi )

    pp = dw.get( 'pfe', self.dwi )                         # External force
    gg = dw.get( 'gfe', self.dwi )
    self.util.integrate( cIdx, cW, pp, gg )
    }

    void BrMPM::MPMMaterial::interpolateParticlesToGrid(MPMDatawarehouseP& dw,
        MPMPatchP& patch) {
    # Interpolate particle mass and momentum to the grid
    cIdx,cW = dw.getMult( ['cIdx','cW'], self.dwi )

    pp = dw.get( 'pm', self.dwi )                          # Mass
    gg = dw.get( 'gm', self.dwi)
    self.util.integrate( cIdx, cW, pp, gg )

    pp = dw.get( 'pw', self.dwi )                          # Momentum
    gg = dw.get( 'gw', self.dwi )
    self.util.integrate( cIdx, cW, pp, gg )
    }

    void BrMPM::MPMMaterial::computeStressTensor(MPMDatawarehouseP& dw,
        MPMPatchP& patch) {
    pf  = dw.get( 'pF', self.dwi )                # Deformation Gradient
    pvs = dw.get( 'pVS', self.dwi )               # Volume * Stress
    pv  = dw.get( 'pVol', self.dwi )              # Volume
    pn  = dw.get( 'pn', self.dwi )              # Volume

    for (ii,pfi,pvi) in izip(count(),pf,pv):
        S,Ja = self.mm.getStress( pfi )     # Get stress and det(pf)
        pn[ii] = Ja
        pvs[ii] = S * pvi * Ja              # Stress * deformed volume
        if not self.ignoreNegJ:
            if Ja < 0:  raise JacobianError('computeStressTensor','Neg J')
    }

    void BrMPM::MPMMaterial::computeInternalForce(MPMDatawarehouseP& dw,
        MPMPatchP& patch) {
    # Compute internal body forces - integrate divergence of stress to grid
    cIdx,cGrad = dw.getMult( ['cIdx','cGrad'], self.dwi )

    pp = dw.get( 'pVS', self.dwi )                          # Stress*Volume
    gg = dw.get( 'gfi', self.dwi)
    self.util.divergence( cIdx, cGrad, pp, gg )
    }

    void BrMPM::MPMMaterial::computeAndIntegrateAcceleration(
        MPMDatawarehouseP& dw, MPMPatchP& patch, double& tol) {
    # Integrate grid acceleration
    dwi = self.dwi
    a_leap = 1. - (patch.it==0) * 0.5             # Initializes leap-frog

    gm = dw.get( 'gm', dwi )                      # Mass
    gw = dw.get( 'gw', dwi )                      # Momentum
    gwc = dw.get('gwc', dwi )
    gfi = dw.get( 'gfi', dwi )                    # Internal Force
    gfe = dw.get( 'gfe', dwi )                    # External Force
    pm = dw.get('pm',dwi)

    cIdx,cW,cGrad = dw.getMult( ['cIdx','cW','cGrad'], self.dwi )
    pfi = dw.get( 'pfi', dwi )
    self.util.interpolate( cIdx, cW, pfi, gfi )

    pfc = dw.get( 'pfc', dwi )
    self.util.interpolate( cIdx, cW, pfc, gfe )
    gv = dw.get( 'gv', dwi )                      # Velocity
    ga = dw.get( 'ga', dwi )

    gm[:] += tol
    gv[:] = (gw+gwc)/gm
    #gv[:] = gw/gm
    ga[:] = a_leap * (gfe+gfi)/gm
    #ga[:] = (gfe+gfi)/gm
    gv[:] += ga*patch.dt
    }

    void BrMPM::MPMMaterial::interpolateToParticlesAndUpdate(
        MPMDatawarehouseP& dw, MPMPatchP& patch) {
    dwi = self.dwi
    cIdx,cW,cGrad = dw.getMult( ['cIdx','cW','cGrad'], self.dwi )

    pvI = dw.get( 'pvI', dwi )
    pxI = dw.get( 'pxI', dwi )
    pGv = dw.get( 'pGv', dwi )
    ga  = dw.get( 'ga', dwi )
    gv  = dw.get( 'gv', dwi )

    self.util.interpolate( cIdx, cW, pvI, ga )
    self.util.interpolate( cIdx, cW, pxI, gv )
    self.util.gradient( cIdx, cGrad, pGv, gv )

    px = dw.get( 'px', dwi )
    pw = dw.get( 'pw', dwi )
    pm = dw.get( 'pm', dwi )
    pF = dw.get( 'pF', dwi )

    #pw += pvI * pm * patch.dt
    pw[:] = pxI * pm
    px[:] += pxI * patch.dt

    self.util.dotAdd( pF, pGv*patch.dt )                # pF += (pGv*dt).pF
    }

