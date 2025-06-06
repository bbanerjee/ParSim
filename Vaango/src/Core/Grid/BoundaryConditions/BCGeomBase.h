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

#ifndef __CORE_GRID_BOUNDARY_CONDITIONS_BCGeomBase_H__
#define __CORE_GRID_BOUNDARY_CONDITIONS_BCGeomBase_H__

#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Point.h>
#include <Core/Grid/BoundaryConditions/BCData.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/BaseIterator.h>
#include <Core/Grid/Variables/Iterator.h>
#include <iterator>
#include <typeinfo>
#include <vector>

namespace Uintah {

/*!

\class BCGeomBase

\ brief Base class for the boundary condition geometry types.

\author John A. Schmidt \n
Department of Mechanical Engineering \n
University of Utah \n
Center for the Simulation of Accidental Fires and Explosions (C-SAFE) \n\n

*/

using std::vector;
using Uintah::IntVector;
using Uintah::Point;

class BCGeomBase
{
public:
  /*  \struct ParticleBoundarySpec
   *  \author Tony Saad
   *  \date   August 2014
   *  \brief  Struct that holds information about particle boundary conditions
   */
  struct ParticleBndSpec
  {

    /*
     *  \brief Particle boundary type wall or inlet
     */
    enum ParticleBndTypeEnum
    {
      WALL,
      INLET,
      NOTSET
    };

    /*
     *  \brief Wall boundary type: Elastic, Inelastic, etc...
     */
    enum ParticleWallTypeEnum
    {
      ELASTIC,
      INELASTIC,
      PARTIALLYELASTIC
    };

    // Default constructor
    ParticleBndSpec()
    {
      ParticleBndSpec(
        ParticleBndSpec::NOTSET, ParticleBndSpec::ELASTIC, 0.0, 0.0);
    }

    // Constructor
    ParticleBndSpec(ParticleBndTypeEnum bndTypet,
                    ParticleWallTypeEnum wallTypet,
                    const double restitution,
                    const double pPerSec)
    {
      bndType         = bndTypet;
      wallType        = wallTypet;
      restitutionCoef = restitution;
      particlesPerSec = pPerSec;
    }

    // Copy constructor
    ParticleBndSpec(const ParticleBndSpec& rhs) = default;
    auto
    operator=(const ParticleBndSpec& rhs) -> ParticleBndSpec& = default;

    /*
     *  \brief Checks whether a particle boundary condition has been specified
     * on this BCGeometry
     */
    [[nodiscard]] auto
    hasParticleBC() const -> bool
    {
      return (bndType != ParticleBndSpec::NOTSET);
    }

    ParticleBndTypeEnum bndType;
    ParticleWallTypeEnum wallType;
    double restitutionCoef;
    double particlesPerSec;
  };

  /// Constructor
  BCGeomBase();

  /// Copy constructor
  BCGeomBase(const BCGeomBase& rhs);

  /// Assignment operator
  auto
  operator=(const BCGeomBase& rhs) -> BCGeomBase&;

  /// Destructor
  virtual ~BCGeomBase();

  /// Equality test
  virtual auto
  operator==(const BCGeomBase&) const -> bool = 0;

  /// Make a clone
  virtual auto
  clone() -> std::shared_ptr<BCGeomBase> = 0;

  /// Get the boundary condition data
  virtual void
  getBCData(BCData& bc) const = 0;

  /// For old boundary conditions
  virtual void
  addBCData(BCData& bc) = 0;

  /// For old boundary conditions
  virtual void
  addBC(BoundCondBaseSP bc) = 0;

  /// Allows a component to add a boundary condition, which already has an
  /// iterator.
  //  This method is exactly the same as addBC, except it applies to two
  //  additional geometries, differences and unions.  Since these geometries are
  //  not "atomic" and can consist of other sub-geometries the function addBC
  //  intentionally will not add boundary conditions to these objects. Hence,
  //  this function forces the boundary conditions to be set regardless of the
  //  inherited object's special properties.
  virtual void
  sudoAddBC(BoundCondBaseSP& bc) = 0;

  void
  getCellFaceIterator(Iterator& b_ptr);

  void
  getNodeFaceIterator(Iterator& b_ptr);

  auto
  hasIterator() -> bool
  {
    return (d_cells.size() > 0);
  }

  /// Determine if a point is inside the geometry where the boundary
  /// condition is applied.
  [[nodiscard]] virtual auto
  inside(const Point& p) const -> bool = 0;

  /// Print out the type of boundary condition -- debugging
  virtual void
  print() = 0;

  /// Determine the cell centered boundary and node centered boundary
  /// iterators.
  virtual void
  determineIteratorLimits(Patch::FaceType face,
                          const Patch* patch,
                          std::vector<Point>& test_pts);

  /*
       \Author  Tony Saad
       \Date    September 2014
       \brif    Determine the iterator associated with this geometry when it is
     used as an interior boundary.
       */
  virtual void
  determineInteriorBndIteratorLimits(const Patch::FaceType face,
                                     const Patch* patch);

  /// Print out the iterators for the boundary.
  void
  printLimits() const;

  /// Get the name for this boundary specification
  auto
  getBCName() -> string
  {
    return d_bcname;
  };
  void
  setBCName(std::string bcname)
  {
    d_bcname = bcname;
  };

  /// Get the type for this boundary specification (type is usually associated
  /// with a user-friendly boundary type such as Wall, Inlet, Outflow...
  auto
  getBndType() -> std::string
  {
    return d_bndtype;
  }
  void
  setBndType(std::string bndType)
  {
    d_bndtype = bndType;
  }

  // Particle-related functionality
  auto
  getParticleBndSpec() -> ParticleBndSpec
  {
    return d_particleBndSpec;
  }

  void
  setParticleBndSpec(const ParticleBndSpec pBndSpec)
  {
    d_particleBndSpec = pBndSpec;
  }

  auto
  hasParticleBC() -> bool
  {
    return d_particleBndSpec.hasParticleBC();
  }

  auto
  surfaceArea() -> double
  {
    return d_surfaceArea;
  }

  auto
  getOrigin() -> Point
  {
    return d_origin;
  }

protected:
  Iterator d_cells;
  Iterator d_nodes;
  std::string d_bcname;
  std::string d_bndtype;
  ParticleBndSpec d_particleBndSpec;
  double d_surfaceArea;
  Point d_origin;
};

template<class T>
class cmp_type
{
public:
  auto
  operator()(const std::shared_ptr<BCGeomBase>& p) -> bool
  {
    return (typeid(T) == typeid(*p.get()));
  }
  auto
  operator()(const BCGeomBase* p) -> bool
  {
    return (typeid(T) == typeid(*p));
  }
};

template<class T>
class not_type
{
public:
  auto
  operator()(const std::shared_ptr<BCGeomBase>& p) -> bool
  {
    return (typeid(T) != typeid(*p.get()));
  }
  auto
  operator()(const BCGeomBase* p) -> bool
  {
    return (typeid(T) != typeid(*p));
  }
};

template<typename T>
class delete_object
{
public:
  void
  operator()(std::shared_ptr<T>& ptr)
  {
    ptr = nullptr;
  }
  void
  operator()(T* ptr)
  {
    delete ptr;
  }
};

} // End namespace Uintah

#endif //__CORE_GRID_BOUNDARY_CONDITIONS_BCGeomBase_H__
