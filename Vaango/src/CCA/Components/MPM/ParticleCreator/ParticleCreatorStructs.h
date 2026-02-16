/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2026 Biswajit Banerjee, Parresia Research Ltd., NZ
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef PARTICLE_CREATOR_STRUCTS_H
#define PARTICLE_CREATOR_STRUCTS_H

#include <Core/Grid/Variables/ParticleVariable.h>

namespace Uintah {

using particleIndex = int; // TODO: Check whether this should be uint64
using particleId    = int;

class GeometryObject;
class Patch;
class DataWarehouse;
class MPMFlags;
class MPMMaterial;
class MPMLabel;
class AMRMPMLabel;
class HydroMPMLabel;
class ParticleSubset;
class VarLabel;

using GeomName   = std::pair<std::string, GeometryObject*>;
using GeomPoint  = std::map<GeometryObject*, std::vector<Point>>;
using GeomScalar = std::map<GeomName, std::vector<double>>;
using GeomVector = std::map<GeomName, std::vector<Vector>>;
using GeomTensor = std::map<GeomName, std::vector<Matrix3>>;
using VecGeometryObjectSP = std::vector<std::shared_ptr<GeometryObject>>;
struct ObjectVars
{
  GeomPoint points;
  GeomScalar scalars;
  GeomVector vectors;
  GeomTensor tensors;
};

// Helper structure to organize particle variable name-object pairs
struct ParticleVarPairs {
  GeomName volume;
  GeomName temperature;
  GeomName color;
  GeomName concentration;
  GeomName poscharge;
  GeomName negcharge;
  GeomName permittivity;
  GeomName externalforce;
  GeomName fiber;
  GeomName velocity;
  GeomName area;
  GeomName size;
  
  explicit ParticleVarPairs(GeometryObject* obj)
    : volume("p.volume", obj)
    , temperature("p.temperature", obj)
    , color("p.color", obj)
    , concentration("p.concentration", obj)
    , poscharge("p.poscharge", obj)
    , negcharge("p.negcharge", obj)
    , permittivity("p.permittivity", obj)
    , externalforce("p.externalforce", obj)
    , fiber("p.fiber", obj)
    , velocity("p.velocity", obj)
    , area("p.area", obj)
    , size("p.size", obj)
  {}
};

struct ParticleVars
{
  ParticleVariable<Point> position;
  ParticleVariable<Point> pX0;
  ParticleVariable<Vector> pDisp, pVelocity, pAcc, pExternalForce;
  ParticleVariable<Matrix3> pSize;
  ParticleVariable<double> pMass, pVolume, pTemperature, pSpecificVolume,
    pErosion;
  ParticleVariable<double> pColor, pTempPrevious, p_q;
  ParticleVariable<long64> pParticleID;
  ParticleVariable<long64> pRigidBodyID;
  ParticleVariable<Vector> pAngularVelocity, pTorque;
  ParticleVariable<Matrix3> pOrientation, pInertiaTensor;
  ParticleVariable<double> pRadius;
  ParticleVariable<Vector> pFiberDir;
  ParticleVariable<int> pLoadCurveID;
  ParticleVariable<IntVector> pLoadCurveIDVector;

  // Body forces
  ParticleVariable<Vector> pBodyForceAcc;
  ParticleVariable<double> pCoriolisImportance;

  // ImplicitParticleCreator
  ParticleVariable<double> pVolumeold;

  // MembraneParticleCreator
  ParticleVariable<Vector> pTang1, pTang2, pNorm;

  // AMR
  ParticleVariable<int> pRefined;
  ParticleVariable<int> pLastLevel;

  // Switch between explicit and implicit MPM
  ParticleVariable<double> pExternalHeatFlux;

  // For friction contact
  ParticleVariable<double> pSurface;
  ParticleVariable<double> pSurfaceGrad;

  // Scalar Diffusion
  ParticleVariable<double> pConcentration;
  ParticleVariable<double> pConcPrevious;
  ParticleVariable<Vector> pConcGrad;
  ParticleVariable<double> pExternalScalarFlux;
  ParticleVariable<double> pPosCharge;
  ParticleVariable<double> pNegCharge;
  ParticleVariable<Vector> pPosChargeGrad;
  ParticleVariable<Vector> pNegChargeGrad;
  ParticleVariable<double> pPermittivity;
  ParticleVariable<Vector> pArea;

  // Hydro-mechanical coupling MPM
  ParticleVariable<double> pFluidMass, pSolidMass, pPorePressure, pPorosity;
  ParticleVariable<Vector> pFluidVelocity, pFluidAcceleration;
  ParticleVariable<Vector> pPrescribedPorePressure;
};

} // end namespace Uintah

#endif // PARTICLE_CREATOR_STRUCTS_H