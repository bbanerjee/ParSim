/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

//  CZMaterial.cc

#include <CCA/Components/MPM/CohesiveZone/CZMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/Core/CZLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Geometry/IntVector.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <iostream>
#include <list>
#include <string>

using namespace Uintah;

// Default constructor
CZMaterial::CZMaterial()
{
  d_mpm_labels = std::make_unique<MPMLabel>();
}

// Standard Constructor
CZMaterial::CZMaterial(ProblemSpecP& ps,
                       const MaterialManagerP& mat_manager,
                       const MPMFlags* flags)
  : Material(ps)
{
  d_mpm_labels = std::make_unique<MPMLabel>();

  // The standard set of initializations needed
  standardInitialization(ps, flags);

  d_cohesive_zone = std::make_unique<CohesiveZone>(
    this, mat_manager, d_mpm_labels.get(), flags);
}

void
CZMaterial::standardInitialization(ProblemSpecP& ps, const MPMFlags* flags)

{
  ps->require("delta_n", d_delta_n);
  ps->require("delta_t", d_delta_t);
  ps->require("sig_max", d_sig_max);
  ps->require("tau_max", d_tau_max);
  ps->require("cz_filename", d_cz_filename);
  ps->getWithDefault("do_rotation", d_do_rotation, false);
  ps->getWithDefault("delta_n_fail", d_delta_n_fail, 4.0 * d_delta_n);
  ps->getWithDefault("delta_t_fail", d_delta_t_fail, 4.0 * d_delta_t);
}

void
CZMaterial::registerParticleState(constVarLabel2DArray& state,
                                  constVarLabel2DArray& state_preReloc)
{
  state.push_back(d_cohesive_zone->returnCohesiveZoneState());
  state_preReloc.push_back(d_cohesive_zone->returnCohesiveZoneStatePreReloc());
}

ProblemSpecP
CZMaterial::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP cz_ps = ps->appendChild("cohesive_zone");

  cz_ps->appendElement("delta_n", d_delta_n);
  cz_ps->appendElement("delta_t", d_delta_t);
  cz_ps->appendElement("sig_max", d_sig_max);
  cz_ps->appendElement("tau_max", d_tau_max);
  cz_ps->appendElement("cz_filename", d_cz_filename);
  cz_ps->appendElement("do_rotation", d_do_rotation);
  cz_ps->appendElement("cz_filename", d_cz_filename);
  cz_ps->appendElement("do_rotation", d_do_rotation);

  return cz_ps;
}

void
CZMaterial::copyWithoutGeom(ProblemSpecP& ps,
                            const CZMaterial* mat,
                            const MPMFlags* flags)
{
  d_delta_n      = mat->d_delta_n;
  d_delta_t      = mat->d_delta_t;
  d_delta_n_fail = mat->d_delta_n_fail;
  d_delta_t_fail = mat->d_delta_t_fail;
  d_sig_max      = mat->d_sig_max;
  d_tau_max      = mat->d_tau_max;
  d_cz_filename  = mat->d_cz_filename;
  d_do_rotation  = mat->d_do_rotation;
}

CohesiveZone*
CZMaterial::getCohesiveZone()
{
  return d_cohesive_zone.get();
}

double
CZMaterial::getCharLengthNormal() const
{
  return d_delta_n;
}

double
CZMaterial::getCharLengthTangential() const
{
  return d_delta_t;
}

double
CZMaterial::getNormalFailureDisplacement() const
{
  return d_delta_n_fail;
}

double
CZMaterial::getTangentialFailureDisplacement() const
{
  return d_delta_t_fail;
}

double
CZMaterial::getCohesiveNormalStrength() const
{
  return d_sig_max;
}

double
CZMaterial::getCohesiveTangentialStrength() const
{
  return d_tau_max;
}

string
CZMaterial::getCohesiveFilename() const
{
  return d_cz_filename;
}

bool
CZMaterial::getDoRotation() const
{
  return d_do_rotation;
}

void
CZMaterial::computeRotationMatrix(Matrix3& Rotation,
                                  Matrix3& Rotation_tang,
                                  const Vector& cznorm,
                                  const Vector czsep) const
{
  double disp  = czsep.length();
  Vector axis  = Cross(cznorm, czsep / disp);
  double theta = std::acos(1.0 - (0.5 * Dot(czsep, czsep)));
  double ca    = std::cos(theta);
  double sa    = std::sin(theta);

  Rotation(0, 0) = (ca - axis[0] * axis[0]) * ca + axis[0] * axis[0];
  Rotation(0, 1) = (-axis[0] * axis[1]) * ca + axis[0] * axis[1] - axis[2] * sa;
  Rotation(0, 2) = (-axis[0] * axis[2]) * ca + axis[0] * axis[2] + axis[1] * sa;
  Rotation(1, 0) = (-axis[1] * axis[0]) * ca + axis[1] * axis[0] + axis[2] * sa;
  Rotation(1, 1) = (ca - axis[1] * axis[1]) * ca + axis[1] * axis[1];
  Rotation(1, 2) = (-axis[1] * axis[2]) * ca + axis[1] * axis[2] - axis[0] * sa;
  Rotation(2, 0) = (-axis[2] * axis[0]) * ca + axis[2] * axis[0] - axis[1] * sa;
  Rotation(2, 1) = (-axis[2] * axis[1]) * ca + axis[2] * axis[1] + axis[0] * sa;
  Rotation(2, 2) = (ca - axis[2] * axis[2]) * ca + axis[2] * axis[2];

  Vector axisx(1.0, 0.0, 0.0);
  Vector axisy(0.0, 1.0, 0.0);
  Vector axisz(0.0, 0.0, 1.0);
  double alpha = 0.;
  double beta  = 0.;
  double gamma = 0.;
  if (std::abs(Rotation(2, 0) != 1)) {
    beta  = -std::asin(Rotation(2, 0));
    alpha = std::atan(Rotation(2, 1) / Rotation(2, 2));
    gamma = std::atan(Rotation(1, 0) / Rotation(0, 0));
  } else {
    gamma = 0;
    alpha = gamma + std::atan(Rotation(0, 1) / Rotation(0, 2));
    beta  = 0.5 * M_PI;
  }
  Matrix3 Rotationx;
  Matrix3 Rotationy;
  Matrix3 Rotationz;
  double calpha = std::cos(alpha);
  double cbeta  = std::cos(beta);
  double salpha = std::sin(alpha);
  double sbeta  = std::sin(beta);
  Rotationx(0, 0) =
    (calpha - axisx[0] * axisx[0]) * calpha + axisx[0] * axisx[0];
  Rotationx(0, 1) =
    (-axisx[0] * axisx[1]) * calpha + axisx[0] * axisx[1] - axisx[2] * salpha;
  Rotationx(0, 2) =
    (-axisx[0] * axisx[2]) * calpha + axisx[0] * axisx[2] + axisx[1] * salpha;
  Rotationx(1, 0) =
    (-axisx[1] * axisx[0]) * calpha + axisx[1] * axisx[0] + axisx[2] * salpha;
  Rotationx(1, 1) =
    (calpha - axisx[1] * axisx[1]) * calpha + axisx[1] * axisx[1];
  Rotationx(1, 2) =
    (-axisx[1] * axisx[2]) * calpha + axisx[1] * axisx[2] - axisx[0] * salpha;
  Rotationx(2, 0) =
    (-axisx[2] * axisx[0]) * calpha + axisx[2] * axisx[0] - axisx[1] * salpha;
  Rotationx(2, 1) =
    (-axisx[2] * axisx[1]) * calpha + axisx[2] * axisx[1] + axisx[0] * salpha;
  Rotationx(2, 2) =
    (calpha - axisx[2] * axisx[2]) * calpha + axisx[2] * axisx[2];
  Rotationy(0, 0) =
    (cbeta - axisy[0] * axisy[0]) * calpha + axisy[0] * axisy[0];
  Rotationy(0, 1) =
    (-axisy[0] * axisy[1]) * cbeta + axisy[0] * axisy[1] - axisy[2] * sbeta;
  Rotationy(0, 2) =
    (-axisy[0] * axisy[2]) * cbeta + axisy[0] * axisy[2] + axisy[1] * sbeta;
  Rotationy(1, 0) =
    (-axisy[1] * axisy[0]) * cbeta + axisy[1] * axisy[0] + axisy[2] * sbeta;
  Rotationy(1, 1) = (cbeta - axisy[1] * axisy[1]) * cbeta + axisy[1] * axisy[1];
  Rotationy(1, 2) =
    (-axisy[1] * axisy[2]) * cbeta + axisy[1] * axisy[2] - axisy[0] * sbeta;
  Rotationy(2, 0) =
    (-axisy[2] * axisy[0]) * cbeta + axisy[2] * axisy[0] - axisy[1] * sbeta;
  Rotationy(2, 1) =
    (-axisy[2] * axisy[1]) * cbeta + axisy[2] * axisy[1] + axisy[0] * sbeta;
  Rotationy(2, 2) = (cbeta - axisy[2] * axisy[2]) * cbeta + axisy[2] * axisy[2];
  Rotationz(0, 0) =
    (cos(gamma) - axisz[0] * axisz[0]) * cos(gamma) + axisz[0] * axisz[0];
  Rotationz(0, 1) = (-axisz[0] * axisz[1]) * cos(gamma) + axisz[0] * axisz[1] -
                    axisz[2] * sin(gamma);
  Rotationz(0, 2) = (-axisz[0] * axisz[2]) * cos(gamma) + axisz[0] * axisz[2] +
                    axisz[1] * sin(gamma);
  Rotationz(1, 0) = (-axisz[1] * axisz[0]) * cos(gamma) + axisz[1] * axisz[0] +
                    axisz[2] * sin(gamma);
  Rotationz(1, 1) =
    (cos(gamma) - axisz[1] * axisz[1]) * cos(gamma) + axisz[1] * axisz[1];
  Rotationz(1, 2) = (-axisz[1] * axisz[2]) * cos(gamma) + axisz[1] * axisz[2] -
                    axisz[0] * sin(gamma);
  Rotationz(2, 0) = (-axisz[2] * axisz[0]) * cos(gamma) + axisz[2] * axisz[0] -
                    axisz[1] * sin(gamma);
  Rotationz(2, 1) = (-axisz[2] * axisz[1]) * cos(gamma) + axisz[2] * axisz[1] +
                    axisz[0] * sin(gamma);
  Rotationz(2, 2) =
    (cos(gamma) - axisz[2] * axisz[2]) * cos(gamma) + axisz[2] * axisz[2];

  Rotation_tang = Rotationz * Rotationy * Rotationx;
}
