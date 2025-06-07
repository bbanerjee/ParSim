/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef __MPM_FLAGS_H__
#define __MPM_FLAGS_H__
#include <CCA/Ports/Output.h>
#include <Core/Grid/MPMInterpolators/ParticleInterpolator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <string>
#include <vector>

namespace Uintah {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class MPMFlags
  \brief A structure that store the flags used for a MPM simulation
  \author Biswajit Banerjee \n
  C-SAFE and Department of Mechanical Engineering \n
  University of Utah \n
*/
/////////////////////////////////////////////////////////////////////////////

class MPMFlags
{

public:
  enum IntegratorType
  {
    Explicit,
    Implicit,
    Fracture
  };

  // flags for turning on/off the reduction variable calculations
  struct ReductionVars
  {
    bool mass;
    bool momentum;
    bool thermalEnergy;
    bool strainEnergy;
    bool accStrainEnergy;
    bool KE;
    bool volDeformed;
    bool centerOfMass;
    bool sumTransmittedForce;
  };

  MPMFlags(const ProcessorGroup* myworld);
  MPMFlags(const MPMFlags& state) = delete;
  MPMFlags&
  operator=(const MPMFlags& state) = delete;
  virtual ~MPMFlags()              = default;

  virtual void
  readMPMFlags(ProblemSpecP& ps, Output* dataArchive);
  virtual void
  outputProblemSpec(ProblemSpecP& ps);
  bool
  doMPMOnLevel(int level, int numLevels) const;

public:
  const ProcessorGroup* d_myworld;
  Output* d_output{ nullptr };

  Uintah::Vector d_gravity;
  IntegratorType d_integrator;
  std::vector<std::string> d_boundaryTractionFaceStrings;
  std::unique_ptr<ReductionVars> d_reductionVars;
  std::unique_ptr<ParticleInterpolator> d_interpolator;

  bool d_withICE{ false };

  bool d_addNewMaterial;             // Flag to decide whether to create
  bool d_axisymmetric;               // Use axisymmetric?
  bool d_artificialViscosity;        // Turn artificial viscosity on/off
  bool d_artificialViscosityHeating; // Include heating due to AV
  bool d_createNewParticles;         // Flag to decide whether to create
                                     // new particles after failure
  bool d_canAddMPMMaterial;          //
  bool d_computeCollinearNormals;    //
  bool d_computeNodalHeatFlux;       // compute the auxilary nodal heat flux
  bool d_computeScaleFactor;         // compute the scale factor for viz
  bool d_deleteRogueParticles;       // Flag to delete rogue particles
  bool d_doContactFriction;          //
  bool d_doErosion;                  // Flag to decide whether to erode or not
  bool d_doExplicitHeatConduction;   //
  bool d_doGridReset;                // Default is true, standard MPM
  bool d_doImplicitHeatConduction;   //
  bool d_doPressureStabilization;    //
  bool d_doScalarDiffusion;          // Flag for scalar diffusion
  bool d_doThermalExpansion;         // Decide whether to do thermExp or not
  bool d_doTransientImplicitHeatConduction; //
  bool d_exactDeformation; // Set steps exactly to match times in prescribed
                           // deformation file
  bool d_fracture;         // to turn on fracture
  bool d_initializeStressFromBodyForce; // Flag for using body force to
                                        // initialize stress
  bool d_insertParticles;               // Activate particles according to color
  bool d_prescribeDeformation; // Prescribe deformation via a table of U and R
  bool d_useCohesiveZones;     // Flag for using cohesive zones
  bool d_useCoordRotation;     // Coordinate rotation on/off
  bool d_useCBDI;       // Flag for using CBDI boundary condition treatment
  bool d_useCPTI;       // Flag for using CPTI interpolation
  bool d_useLoadCurves; // Flag for using load curves
  bool d_useLoadCurvesVector{ false }; // Flag for using load curves
  bool d_withColor;                    // to turn on the color variable
  bool d_AMR;                          // Do AMR?
  bool d_GEVelProj;                    // Flag for adaptive mesh refinement
  bool d_refineParticles;              // Flag for refinement
  bool d_useXPIC{ false };             // Use XPIC (Nairn et al.) algorithm
  bool d_updateStressLast{ false };
  bool d_deleteGeometryObjects{ false };

  // For hydro-mechanical coupling
  bool d_coupledFlow{ false };

  // For scalar diffusion
  double d_autoCycleMin{ 0.1 };
  double d_autoCycleMax{ 0.9 };
  bool d_doAutoCycleBC{ false };
  bool d_autoCycleUseMinMax{ false };

  int d_ndim{ 3 };
  int d_8or27;        // Number of nodes a particle can interact with
  int d_minGridLevel; // Only do MPM on this grid level
  int d_maxGridLevel; // Only do MPM on this grid level
  int d_numTermsSeriesDefGrad{ 1 }; // Number of terms in series expansion
                                    // for deformation gradient calculation
  int d_extraSolverFlushes;         // Have PETSc do more flushes to save memory

  double d_addFrictionWork;      // 1 == add , 0 == do not add
  double d_artificialDampCoeff;  //
  double d_artificialViscCoeff1; // Artificial viscosity coefficient 1
  double d_artificialViscCoeff2; // Artificial viscosity coefficient 2
  double d_cpdiLcrit;            // for cpdi interpolator maximum fractional
                                 // cell size for a particle
  double d_forceIncrementFactor; //

  double d_maxVel;      // Maxmimum particle velocity before  deletion
  double d_minPartMass; // Minimum particle mass before deletion
  double d_minMassForAcceleration; // Minimum mass to allow division by in
                                   // computing acceleration

  std::string d_interpolatorType; // Type of particle-grid interaction
  std::string d_integratorType;   // Explicit or implicit time integration
  std::string d_defGradAlgorithm; // Deformation gradient algorithm
  std::string d_erosionAlgorithm; // Algorithm to erode material points
  std::string d_prescribedDeformationFile; //
  std::string d_insertParticlesFile;       // File containing activation plan
  std::string d_mmsType;                   // MMS Flag

  // For rotating coordinate system
  double d_coordRotationSpeed;               // Speed of rotation
  Uintah::Point d_coordRotationCenter;       // Center of rotation
  Uintah::Point d_coordRotationBodyRefPoint; // Reference point in rotating body
  Uintah::Vector d_coordRotationAxis;        // Axis of rotation

  // So all components can know how many particle ghost cells to ask for
  Ghost::GhostType particle_ghost_type{ Ghost::None };
  int particle_ghost_layer{ 0 };
};

} // End namespace Uintah

#endif // __MPM_FLAGS_H__
