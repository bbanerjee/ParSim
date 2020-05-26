/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/MPMFlags.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/AxiGIMPInterpolator.h>
#include <Core/Grid/AxiLinearInterpolator.h>
#include <Core/Grid/GIMPInterpolator.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/LinearInterpolator.h>
#include <Core/Grid/axiCpdiInterpolator.h>
#include <Core/Grid/cpdiInterpolator.h>
#include <Core/Grid/cptiInterpolator.h>
//#include <Core/Grid/fastCpdiInterpolator.h>
//#include <Core/Grid/fastAxiCpdiInterpolator.h>
#include <Core/Grid/BSplineInterpolator.h>
#include <Core/Grid/TOBSplineInterpolator.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Util/DebugStream.h>
#include <iostream>

using namespace Uintah;

static DebugStream dbg("MPMFlags", false);

MPMFlags::MPMFlags(const ProcessorGroup* myworld)
{
  d_myworld = myworld;
  d_withICE = false;

  d_gravity = Vector(0., 0., 0.);

  d_reductionVars                  = std::make_unique<ReductionVars>();
  d_reductionVars->mass            = false;
  d_reductionVars->momentum        = false;
  d_reductionVars->thermalEnergy   = false;
  d_reductionVars->strainEnergy    = false;
  d_reductionVars->accStrainEnergy = false;
  d_reductionVars->KE              = false;
  d_reductionVars->volDeformed     = false;
  d_reductionVars->centerOfMass    = false;

  d_integratorType = "explicit";
  d_integrator     = Explicit;

  d_interpolator     = std::make_unique<LinearInterpolator>();
  d_interpolatorType = "linear";
  d_useCBDI          = false;
  d_useCPTI          = false;

  d_axisymmetric = false;
  d_doGridReset = true;
  d_withICE    = false;

  d_minGridLevel = 0;
  d_maxGridLevel = 1000;

  d_minPartMass             = 3.e-15;
  d_minMassForAcceleration = 0;
  d_maxVel                   = 3.e105;

  d_artificialViscosity        = false;
  d_artificialViscosityHeating = false;
  d_artificialViscCoeff1       = 0.2;
  d_artificialViscCoeff2       = 2.0;
  d_artificialDampCoeff        = 0.0;

  d_useLoadCurves    = false;
  d_useCohesiveZones = false;

  d_insertParticles    = false;
  d_createNewParticles = false;
  d_addNewMaterial     = false;
  d_withColor          = false;
  d_canAddMPMMaterial  = false;

  d_erosionAlgorithm     = "none";
  d_doErosion            = false;
  d_deleteRogueParticles = false;
  d_fracture             = false;

  d_forceIncrementFactor    = 1.0;
  d_doPressureStabilization = false;

  d_doThermalExpansion      = true;
  d_doContactFriction       = false;
  d_addFrictionWork         = 0.0; // don't do frictional heating by default
  d_computeCollinearNormals = false;

  d_extraSolverFlushes       = 0; // Have PETSc do more flushes to save memory
  d_doImplicitHeatConduction = false;
  d_doTransientImplicitHeatConduction = true;
  d_doExplicitHeatConduction          = true;
  d_computeNodalHeatFlux              = false;

  d_computeScaleFactor        = false;
  d_prescribeDeformation      = false;
  d_prescribedDeformationFile = "time_defgrad_rotation";
  d_exactDeformation          = false;

  // MMS
  if (d_mmsType == "AxisAligned") {
    d_mmsType = "AxisAligned";
  } else if (d_mmsType == "GeneralizedVortex") {
    d_mmsType = "GeneralizedVortex";
  } else if (d_mmsType == "ExpandingRing") {
    d_mmsType = "ExpandingRing";
  } else if (d_mmsType == "AxisAligned3L") {
    d_mmsType = "AxisAligned3L";
  }

  // Deformation gradient computer
  d_defGradAlgorithm      = "first_order";
  d_numTermsSeriesDefGrad = 1;

  // Rotating coordinate system
  d_useCoordRotation    = false;
  d_coordRotationCenter = Uintah::Point(0.0, 0.0, 0.0);  // Default is origin
  d_coordRotationAxis   = Uintah::Vector(1.0, 0.0, 0.0); // Default is x-axis
  d_coordRotationSpeed  = 0.0; // **NOTE** Scalar angular speed
  d_coordRotationBodyRefPoint = Uintah::Point(0.0, 0.0, 0.0); // Reference point
  // in rotating body

  // Initialize stress using body force (lithostatic)
  d_initializeStressFromBodyForce = false;

  // Initialize scalar diffusion flag
  d_doScalarDiffusion = false;

  // Initialize AMR flags
  d_AMR             = false;
  d_GEVelProj       = false;
  d_refineParticles = false;
}

void
MPMFlags::readMPMFlags(ProblemSpecP& ps, Output* dataArchive)
{
  ProblemSpecP root         = ps->getRootNode();
  ProblemSpecP mpm_flag_ps  = root->findBlock("MPM");
  ProblemSpecP phys_cons_ps = root->findBlock("PhysicalConstants");

  if (phys_cons_ps) {
    phys_cons_ps->require("gravity", d_gravity);
  } else if (mpm_flag_ps) {
    mpm_flag_ps->require("gravity", d_gravity);
  } else {
    d_gravity = Vector(0, 0, 0);
  }

  //__________________________________
  //  Set the on/off flags to determine which
  // reduction variables are computed
  d_reductionVars->mass          = dataArchive->isLabelSaved("TotalMass");
  d_reductionVars->momentum      = dataArchive->isLabelSaved("TotalMomentum");
  d_reductionVars->thermalEnergy = dataArchive->isLabelSaved("ThermalEnergy");
  d_reductionVars->KE            = dataArchive->isLabelSaved("KineticEnergy");
  d_reductionVars->strainEnergy  = dataArchive->isLabelSaved("StrainEnergy");
  d_reductionVars->accStrainEnergy =
    dataArchive->isLabelSaved("AccStrainEnergy");
  d_reductionVars->volDeformed =
    dataArchive->isLabelSaved("TotalVolumeDeformed");
  d_reductionVars->centerOfMass =
    dataArchive->isLabelSaved("CenterOfMassPosition");

  if (!mpm_flag_ps)
    return;

  mpm_flag_ps->get("boundary_traction_faces", d_bndyFaceTxtList);

  mpm_flag_ps->get("time_integrator", d_integratorType);
  if (d_integratorType == "implicit") {
    d_integrator = Implicit;
  } else if (d_integratorType == "fracture") {
    d_integrator = Fracture;
    d_fracture   = true;
  } else {
    d_integrator = Explicit;
  }
  int junk = 0;
  mpm_flag_ps->get("nodes8or27", junk);
  if (junk != 0) {
    std::cerr << "nodes8or27 is deprecated, use " << "\n";
    std::cerr << "<interpolator>type</interpolator>" << "\n";
    std::cerr << "where type is one of the following:" << "\n";
    std::cerr << "linear, gimp, 3rdorderBS, cpdi, cpti" << "\n";
    exit(1);
  }

  mpm_flag_ps->get("interpolator", d_interpolatorType);
  mpm_flag_ps->getWithDefault("cpdi_lcrit", d_cpdiLcrit, 1.e10);
  mpm_flag_ps->get("axisymmetric", d_axisymmetric);
  mpm_flag_ps->get("with_color", d_withColor);
  mpm_flag_ps->get("artificial_damping_coeff", d_artificialDampCoeff);
  mpm_flag_ps->get("artificial_viscosity", d_artificialViscosity);
  if (d_artificialViscosity) {
    d_artificialViscosityHeating = true;
  }
  mpm_flag_ps->get("artificial_viscosity_heating",
                   d_artificialViscosityHeating);
  mpm_flag_ps->get("artificial_viscosity_coeff1", d_artificialViscCoeff1);
  mpm_flag_ps->get("artificial_viscosity_coeff2", d_artificialViscCoeff2);
  mpm_flag_ps->get("use_load_curves", d_useLoadCurves);
  mpm_flag_ps->get("use_CBDI_boundary_condition", d_useCBDI);
  mpm_flag_ps->get("exact_deformation", d_exactDeformation);
  mpm_flag_ps->get("use_cohesive_zones", d_useCohesiveZones);

  if (d_artificialViscosity && d_integratorType == "implicit") {
    if (d_myworld->myrank() == 0) {
      std::cerr << "artificial viscosity is not implemented" << "\n";
      std::cerr << "with implicit time integration" << "\n";
    }
  }

  if (!d_artificialViscosity && d_artificialViscosityHeating) {
    ostringstream warn;
    warn << "ERROR:MPM: You can't have heating due to artificial viscosity "
         << "if artificial_viscosity is not enabled."
         << "\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  mpm_flag_ps->get("forceBC_force_increment_factor", d_forceIncrementFactor);
  mpm_flag_ps->get("create_new_particles", d_createNewParticles);
  mpm_flag_ps->get("manual_new_material", d_addNewMaterial);
  mpm_flag_ps->get("can_add_MPM_material", d_canAddMPMMaterial);
  mpm_flag_ps->get("do_implicit_heat_conduction", d_doImplicitHeatConduction);
  mpm_flag_ps->get("do_transient_implicit_heat_conduction",
                   d_doTransientImplicitHeatConduction);
  mpm_flag_ps->get("do_explicit_heat_conduction", d_doExplicitHeatConduction);
  mpm_flag_ps->get("do_pressure_stabilization", d_doPressureStabilization);
  mpm_flag_ps->get("do_thermal_expansion", d_doThermalExpansion);
  mpm_flag_ps->get("do_grid_reset", d_doGridReset);
  mpm_flag_ps->get("minimum_particle_mass", d_minPartMass);
  mpm_flag_ps->getWithDefault(
    "minimum_mass_for_acc", d_minMassForAcceleration, 1.0e-199);
  mpm_flag_ps->get("maximum_particle_velocity", d_maxVel);
  mpm_flag_ps->get("use_prescribed_deformation", d_prescribeDeformation);
  if (d_prescribeDeformation) {
    mpm_flag_ps->get("prescribed_deformation_file",
                     d_prescribedDeformationFile);
  }
  // MMS
  mpm_flag_ps->get("run_MMS_problem", d_mmsType);
  // Flag for CPTI interpolator
  if (d_interpolatorType == "cpti") {
    d_useCPTI = true;
  }

  mpm_flag_ps->get("insert_particles", d_insertParticles);
  if (d_insertParticles) {
    mpm_flag_ps->require("insert_particles_file", d_insertParticlesFile);
  }

  mpm_flag_ps->get("do_contact_friction_heating", d_doContactFriction);
  if (!d_doContactFriction)
    d_addFrictionWork = 0.0;

  mpm_flag_ps->getWithDefault(
    "collinear_bimaterial_contact_normals", d_computeCollinearNormals, false);

  ProblemSpecP erosion_ps = mpm_flag_ps->findBlock("erosion");
  if (erosion_ps) {
    if (erosion_ps->getAttribute("algorithm", d_erosionAlgorithm)) {
      if (d_erosionAlgorithm == "none")
        d_doErosion = false;
      else
        d_doErosion = true;
    }
  }

  mpm_flag_ps->get("delete_rogue_particles", d_deleteRogueParticles);

  ProblemSpecP DA_ps = root->findBlock("DataArchiver");
  if (DA_ps) {
    for (ProblemSpecP label_iter = DA_ps->findBlock("save"); label_iter != 0;
         label_iter              = label_iter->findNextBlock("save")) {
      map<string, string> labelName;

      label_iter->getAttributes(labelName);
      if (labelName["label"] == "g.heatflux") {
        d_computeNodalHeatFlux = true;
      }
      if (labelName["label"] == "p.scalefactor") {
        d_computeScaleFactor = true;
      }
    }
  }

  ProblemSpecP da_ps = root->findBlock("DataAnalysis");

  if (da_ps) {
    for (ProblemSpecP module_ps = da_ps->findBlock("Module"); module_ps != 0;
         module_ps              = module_ps->findNextBlock("Module")) {
      if (module_ps) {
        map<string, string> attributes;
        module_ps->getAttributes(attributes);
        if (attributes["name"] == "flatPlate_heatFlux") {
          d_computeNodalHeatFlux = true;
        }
      }
    }
  }
  // restart problem spec
  mpm_flag_ps->get("compute_nodal_heat_flux", d_computeNodalHeatFlux);
  mpm_flag_ps->get("compute_scale_factor", d_computeScaleFactor);

  if (d_interpolatorType == "linear") {
    if (d_axisymmetric) {
      d_interpolator = std::make_unique<AxiLinearInterpolator>();
    } else {
      d_interpolator = std::make_unique<LinearInterpolator>();
    }
  } else if (d_interpolatorType == "gimp") {
    if (d_axisymmetric) {
      d_interpolator = std::make_unique<AxiGIMPInterpolator>();
    } else {
      d_interpolator = std::make_unique<GIMPInterpolator>();
    }
  } else if (d_interpolatorType == "3rdorderBS") {
    if (!d_axisymmetric) {
      d_interpolator = std::make_unique<TOBSplineInterpolator>();
    } else {
      ostringstream warn;
      warn << "ERROR:MPM: invalid interpolation type (" << d_interpolatorType
           << ") Can't be used with axisymmetry at this time \n\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }
  } else if (d_interpolatorType == "4thorderBS") {
    if (!d_axisymmetric) {
      d_interpolator = std::make_unique<BSplineInterpolator>();
    } else {
      ostringstream warn;
      warn << "ERROR:MPM: invalid interpolation type (" << d_interpolatorType
           << ") Can't be used with axisymmetry at this time \n\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }
  } else if (d_interpolatorType == "cpdi") {
    if (d_axisymmetric) {
      d_interpolator = std::make_unique<axiCpdiInterpolator>();
    } else {
      d_interpolator = std::make_unique<cpdiInterpolator>();
      d_interpolator->setLcrit(d_cpdiLcrit);
    }
  } else if (d_interpolatorType == "cpti") {
    if (!d_axisymmetric) {
      d_interpolator = std::make_unique<cptiInterpolator>();
      d_interpolator->setLcrit(d_cpdiLcrit);
    }
  }
#if 0
 else if(d_interpolatorType=="fastcpdi"){
    if(d_axisymmetric){
      d_interpolator = std::make_unique<fastAxiCpdiInterpolator>();
    } else{
      d_interpolator = std::make_unique<fastCpdiInterpolator>();
    }
  }
#endif
  else {
    ostringstream warn;
    warn << "ERROR:MPM: invalid interpolation type (" << d_interpolatorType
         << ")"
         << "Valid options are: \n"
         << "linear\n"
         << "gimp\n"
         << "cpdi\n"
         << "cpti\n"
         << "3rdorderBS\n"
         << "4thorderBS\n"
         << "\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  // Get the size of the vectors associated with the interpolator
  d_8or27 = d_interpolator->size();

  mpm_flag_ps->get("extra_solver_flushes", d_extraSolverFlushes);


  // Deformation gradient computer options
  // Options other than "first_order"/"subcycling" should be defined in
  // DeformationGradientComputer
  // e.g. "constant_velgrad", "linear_velgrad", "constant_velgrad_subcycling",
  // etc.
  ProblemSpecP defgrad_ps = mpm_flag_ps->findBlock("deformation_gradient");
  if (defgrad_ps) {
    if (defgrad_ps->getAttribute("algorithm", d_defGradAlgorithm)) {
      if (d_defGradAlgorithm == "first_order") {
        d_numTermsSeriesDefGrad = 1;
      } else if (d_defGradAlgorithm == "subcycling") {
        d_numTermsSeriesDefGrad = 1;
      } else {
        defgrad_ps->get("num_terms", d_numTermsSeriesDefGrad);
      }
    }
  }

  // For rotating coordinate system
  ProblemSpecP coordRotation_ps =
    mpm_flag_ps->findBlock("rotating_coordinate_system");
  if (coordRotation_ps) {
    d_useCoordRotation = true;
    coordRotation_ps->get("rotation_center", d_coordRotationCenter);
    coordRotation_ps->get("rotation_axis", d_coordRotationAxis);
    coordRotation_ps->get("rotation_speed_angular", d_coordRotationSpeed);
    coordRotation_ps->get("body_reference_point", d_coordRotationBodyRefPoint);

    // Do checks
    if (!(d_coordRotationAxis.length2() > 0.0)) {
      ostringstream warn;
      warn << "ERROR:MPM: Rotation axis has zero length: "
           << d_coordRotationAxis << "\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    } else {
      d_coordRotationAxis.normalize(); // Make sure that the rotation axis has
                                       // unit length
    }

    if (d_coordRotationSpeed < 0.0) {
      ostringstream warn;
      warn << "ERROR:MPM: Rotation speed " << d_coordRotationSpeed << " is < 0"
           << "\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }
  }

  // Initialize stress and deformation gradient using body force ?
  mpm_flag_ps->getWithDefault("initialize_stress_using_body_force",
                              d_initializeStressFromBodyForce,
                              false);

  // Scalar diffusion
  mpm_flag_ps->get("do_scalar_diffusion", d_doScalarDiffusion);

  // AMR
  mpm_flag_ps->get("AMR", d_AMR);
  mpm_flag_ps->getWithDefault(
    "use_gradient_enhanced_velocity_projection", d_GEVelProj, false);
  mpm_flag_ps->get("refine_particles", d_refineParticles);

  if (dbg.active()) {
    dbg << "---------------------------------------------------------\n";
    dbg << "MPM Flags "
        << "\n";
    dbg << "---------------------------------------------------------\n";
    dbg << " Time Integration            = " << d_integratorType << "\n";
    dbg << " Interpolation type          = " << d_interpolatorType << "\n";
    dbg << " With Color                  = " << d_withColor << "\n";
    dbg << " Artificial Damping Coeff    = " << d_artificialDampCoeff << "\n";
    dbg << " Artificial Viscosity On     = " << d_artificialViscosity << "\n";
    dbg << " Artificial Viscosity Htng   = " << d_artificialViscosityHeating
        << "\n";
    dbg << " Artificial Viscosity Coeff1 = " << d_artificialViscCoeff1 << "\n";
    dbg << " Artificial Viscosity Coeff2 = " << d_artificialViscCoeff2 << "\n";
    dbg << " Create New Particles        = " << d_createNewParticles << "\n";
    dbg << " Add New Material            = " << d_addNewMaterial << "\n";
    dbg << " RefineParticles             = " << d_refineParticles << "\n";
    dbg << " Delete Rogue Particles?     = " << d_deleteRogueParticles << "\n";
    dbg << " Use Load Curves             = " << d_useLoadCurves << "\n";
    dbg << " Use CBDI boundary condition = " << d_useCBDI << "\n";
    dbg << " Use Cohesive Zones          = " << d_useCohesiveZones << "\n";
    dbg << " ForceBC increment factor    = " << d_forceIncrementFactor << "\n";
    dbg << " Contact Friction Heating    = " << d_addFrictionWork << "\n";
    dbg << " Extra Solver flushes        = " << d_extraSolverFlushes << "\n";
    dbg << "---------------------------------------------------------\n";
  }
}

void
MPMFlags::outputProblemSpec(ProblemSpecP& ps)
{

  ps->appendElement("gravity", d_gravity);
  ps->appendElement("boundary_traction_faces", d_bndyFaceTxtList);

  ps->appendElement("manual_new_material", d_addNewMaterial);
  ps->appendElement("axisymmetric", d_axisymmetric);
  ps->appendElement("artificial_viscosity", d_artificialViscosity);
  ps->appendElement("artificial_viscosity_heating",
                    d_artificialViscosityHeating);
  ps->appendElement("create_new_particles", d_createNewParticles);
  ps->appendElement("can_add_MPM_material", d_canAddMPMMaterial);
  ps->appendElement("collinear_bimaterial_contact_normals",
                    d_computeCollinearNormals);
  ps->appendElement("compute_nodal_heat_flux", d_computeNodalHeatFlux);
  ps->appendElement("compute_scale_factor", d_computeScaleFactor);
  ps->appendElement("delete_rogue_particles", d_deleteRogueParticles);
  ps->appendElement("do_contact_friction_heating", d_doContactFriction);
  ps->appendElement("do_explicit_heat_conduction", d_doExplicitHeatConduction);
  ps->appendElement("do_grid_reset", d_doGridReset);
  ps->appendElement("do_implicit_heat_conduction", d_doImplicitHeatConduction);
  ps->appendElement("do_pressure_stabilization", d_doPressureStabilization);
  ps->appendElement("do_scalar_diffusion", d_doScalarDiffusion);
  ps->appendElement("do_thermal_expansion", d_doThermalExpansion);
  ps->appendElement("do_transient_implicit_heat_conduction",
                    d_doTransientImplicitHeatConduction);
  ps->appendElement("exact_deformation", d_exactDeformation);
  ps->appendElement("initialize_stress_using_body_force",
                    d_initializeStressFromBodyForce);
  ps->appendElement("insert_particles", d_insertParticles);
  ps->appendElement("use_prescribed_deformation", d_prescribeDeformation);
  ps->appendElement("use_cohesive_zones", d_useCohesiveZones);
  ps->appendElement("use_CBDI_boundary_condition", d_useCBDI);
  ps->appendElement("use_load_curves", d_useLoadCurves);
  ps->appendElement("with_color", d_withColor);
  ps->appendElement("AMR", d_AMR);
  ps->appendElement("use_gradient_enhanced_velocity_projection", d_GEVelProj);
  ps->appendElement("refine_particles", d_refineParticles);

  // PETSc
  ps->appendElement("extra_solver_flushes", d_extraSolverFlushes);

  ps->appendElement("artificial_damping_coeff", d_artificialDampCoeff);
  ps->appendElement("artificial_viscosity_coeff1", d_artificialViscCoeff1);
  ps->appendElement("artificial_viscosity_coeff2", d_artificialViscCoeff2);
  ps->appendElement("cpdi_lcrit", d_cpdiLcrit);
  ps->appendElement("forceBC_force_increment_factor", d_forceIncrementFactor);
  ps->appendElement("maximum_particle_velocity", d_maxVel);
  ps->appendElement("minimum_particle_mass", d_minPartMass);
  ps->appendElement("minimum_mass_for_acc", d_minMassForAcceleration);

  ps->appendElement("interpolator", d_interpolatorType);
  ps->appendElement("time_integrator", d_integratorType);

  ProblemSpecP defgrad_ps = ps->appendChild("deformation_gradient");
  defgrad_ps->setAttribute("algorithm", d_defGradAlgorithm);
  defgrad_ps->appendElement("num_terms", d_numTermsSeriesDefGrad);

  ProblemSpecP erosion_ps = ps->appendChild("erosion");
  erosion_ps->setAttribute("algorithm", d_erosionAlgorithm);

  if (d_prescribeDeformation) {
    ps->appendElement("prescribed_deformation_file",
                      d_prescribedDeformationFile);
  }
  if (d_insertParticles) {
    ps->appendElement("insert_particles_file", d_insertParticlesFile);
  }

  ps->appendElement("run_MMS_problem", d_mmsType);

  // Rotating coordinate system
  ProblemSpecP coordRotation_ps = ps->appendChild("rotating_coordinate_system");
  coordRotation_ps->appendElement("rotation_center", d_coordRotationCenter);
  coordRotation_ps->appendElement("rotation_axis", d_coordRotationAxis);
  coordRotation_ps->appendElement("rotation_speed_angular",
                                  d_coordRotationSpeed);
  coordRotation_ps->appendElement("body_reference_point",
                                  d_coordRotationBodyRefPoint);
}

bool
MPMFlags::doMPMOnLevel(int level, int numLevels) const
{
  return (level >= d_minGridLevel && level <= d_maxGridLevel) ||
         (d_minGridLevel < 0 && level == numLevels + d_minGridLevel);
}
