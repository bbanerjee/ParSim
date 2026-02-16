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

#include <Core/GeometryPiece/SpecialGeomPiece.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <format>

namespace Uintah {

const std::string SpecialGeomPiece::TYPE_NAME = "special_geom";

//////////////////////////////////////////////////////////////////////
/* Returns the vector containing the set of particle locations */
//////////////////////////////////////////////////////////////////////
double
SpecialGeomPiece::volume() const
{
  double total_vol = 0;
  if (d_scalars.count("p.volume")) {
    for (auto vol : d_scalars.at("p.volume")) {
      total_vol += vol;
    }
  } else if (d_tensors.count("p.size")) {
    for (auto size : d_tensors.at("p.size")) {
      total_vol += std::abs(size.Determinant());
    }
  }
  return total_vol;
}

Point
SpecialGeomPiece::getCenter() const
{
  Vector center(0, 0, 0);
  double total_vol = 0;
  if (d_scalars.count("p.volume")) {
    const auto& vols = d_scalars.at("p.volume");
    for (size_t i = 0; i < d_points.size(); ++i) {
      center    = center + d_points[i].asVector() * vols[i];
      total_vol += vols[i];
    }
  } else if (d_tensors.count("p.size")) {
    const auto& sizes = d_tensors.at("p.size");
    for (size_t i = 0; i < d_points.size(); ++i) {
      double vol = std::abs(sizes[i].Determinant());
      center    = center + d_points[i].asVector() * vol;
      total_vol += vol;
    }
  } else {
    for (const auto& point : d_points) {
      center = center + point.asVector();
    }
    total_vol = static_cast<double>(d_points.size());
  }

  if (total_vol > 1e-12) {
    center = center / total_vol;
  }
  return Point(center);
}

Matrix3
SpecialGeomPiece::getInertiaTensor() const
{
  Matrix3 I(0.0);
  Point center = getCenter();
  Matrix3 identity;
  identity.Identity();

  if (d_scalars.count("p.volume")) {
    const auto& vols = d_scalars.at("p.volume");
    for (size_t i = 0; i < d_points.size(); ++i) {
      Vector d       = d_points[i] - center;
      double d_dot_d = d.length2();
      I += (identity * d_dot_d - Matrix3(d, d)) * vols[i];
    }
  } else if (d_tensors.count("p.size")) {
    const auto& sizes = d_tensors.at("p.size");
    for (size_t i = 0; i < d_points.size(); ++i) {
      double vol     = std::abs(sizes[i].Determinant());
      Vector d       = d_points[i] - center;
      double d_dot_d = d.length2();
      I += (identity * d_dot_d - Matrix3(d, d)) * vol;
    }
  } else {
    for (const auto& point : d_points) {
      Vector d       = point - center;
      double d_dot_d = d.length2();
      I += (identity * d_dot_d - Matrix3(d, d)); // Unit volume
    }
  }
  return I;
}

std::vector<Point>*
SpecialGeomPiece::getPoints()
{
  return &d_points;
}

//////////////////////////////////////////////////////////////////////
/*! Returns the vector containing the set of particle data using name */
//////////////////////////////////////////////////////////////////////
const SpecialGeomPiece::ParticleScalarData*
SpecialGeomPiece::getScalar(const std::string& scalar_name) const
{
  auto it = d_scalars.find(scalar_name);
  if (it != d_scalars.end()) {
    return &(it->second);
  } else {
    return nullptr;
  }
}

const SpecialGeomPiece::ParticleVectorData*
SpecialGeomPiece::getVector(const std::string& vector_name) const
{
  auto it = d_vectors.find(vector_name);
  if (it != d_vectors.end()) {
    return &(it->second);
  } else {
    return nullptr;
  }
}

const SpecialGeomPiece::ParticleTensorData*
SpecialGeomPiece::getTensor(const std::string& tensor_name) const
{
  auto it = d_tensors.find(tensor_name);
  if (it != d_tensors.end()) {
    return &(it->second);
  } else {
    return nullptr;
  }
}

//////////////////////////////////////////////////////////////////////
/* Deletes the vector containing the set of particle locations */
//////////////////////////////////////////////////////////////////////
void
SpecialGeomPiece::deletePoints()
{
  d_points.clear();
}

//////////////////////////////////////////////////////////////////////
/* Deletes the vector containing the set of particle scalars */
//////////////////////////////////////////////////////////////////////
void
SpecialGeomPiece::deleteScalar(const std::string& name)
{
  try {
    d_scalars.at(name).clear();
  } catch (std::out_of_range& err) {
  }
}

//////////////////////////////////////////////////////////////////////
/* Deletes the vector containing the set of particle vectors          */
//////////////////////////////////////////////////////////////////////
void
SpecialGeomPiece::deleteVector(const std::string& name)
{
  try {
    d_vectors.at(name).clear();
  } catch (std::out_of_range& err) {
  }
}

//////////////////////////////////////////////////////////////////////
/* Deletes the vector containing the set of particle tensors          */
//////////////////////////////////////////////////////////////////////
void
SpecialGeomPiece::deleteTensor(const std::string& name)
{
  try {
    d_tensors.at(name).clear();
  } catch (std::out_of_range& err) {
  }
}

void
SpecialGeomPiece::writePoints(const std::string& f_name, const std::string& var)
{
  if (var == "p.volume") {
    std::ofstream file(f_name.c_str());
    file.setf(std::ios::scientific, std::ios::floatfield);
    file.precision(8);
    file << "x_coord "
         << " y_coord "
         << " z_coord " << var << "\n";
    for (unsigned int jj = 0; jj < d_points.size(); ++jj) {
      file << d_points[jj].x() << " " << d_points[jj].y() << " "
           << d_points[jj].z() << " " << (d_scalars.at(var))[jj] << "\n";
    }
    file.close();
    std::cout << "Wrote output file " << f_name << std::endl;
  } else {
    std::cout << "** No output file " << f_name << " written." << std::endl;
  }
}

int
SpecialGeomPiece::returnPointCount() const
{
  return d_points.size();
}

void
SpecialGeomPiece::setParticleSpacing(double dx)
{
  d_dx = dx;
}

void
SpecialGeomPiece::setCellSize(Vector DX)
{
  d_DX = DX;
}

} // end namespace Uintah