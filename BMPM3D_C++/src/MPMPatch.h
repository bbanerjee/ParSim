#ifndef __MATITI_MPMPATCH_H__
#define __MATITI_MPMPATCH_H__

#include <MPMPatch.h>
#include <Types.h>
#include <MPMDataTypes.h>

//#include <VelocityBCSPArray.h>

#include <GeometryMath/Point3D.h>
#include <GeometryMath/Vector3D.h>
#include <GeometryMath/IntVector3D.h>

#include <ShapeFunctions/MPMShapeFunctionP.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <iostream>

namespace BrMPM {

class MPMPatch
{

public:

  //    friend std::ostream& operator<<(std::ostream& out, const BrMPM::MPMPatch& patch);

public:

  MPMPatch() ;
  ~MPMPatch();

  void initialize(const Uintah::ProblemSpecP& ps);

  const Vector3D& nGhost() {return d_num_ghost;}
  const IntVector3D& nC()  {return d_node_counts;}
  const IntVector3D& ppe() {return d_num_particles_per_cell;}
  const Vector3D& dX()  {return d_cell_size;}
  const Point3D& x0() {return d_lower;}
  double dt() {return d_delT;}
  int it() {return d_iteration;}
  double getTolerance() {return d_tol;}

  // Shape functions have to be created in the initialize stage
  const MPMShapeFunctionP& shape() {return d_shape;}

  // Initialize grid
  void initGrid(DoubleNodeData& gx);

  bool insidePatch(const Point3D& point) const;
  bool allInsidePatch(const Point3DParticleData& points) const;

  inline void stepTime() {d_time += d_delT; d_iteration += 1;}


  //const Vector3D& ghost() const {return d_num_ghost;}
  //   const double& thick() const {return d_thick;}
  //const int& particlesperelement() const {return d_num_particles_per_cell;}
  // const IntArray3& numGrids() const;
  //const double totalGrids() const;

  //const Point3D lower() const {return d_lower;}
  //const Point3D upper() const {return d_upper;}

  //const Vector3D& cellSize()  {return d_cell_size;}

  //const Vector3D& numGrids()  {return d_num_grids;}
  //const std::vector<Point3D> gridsPosition()  {return d_gridsPosition;}

  //void findGradeIndex(const Point3D& point,
  //    IntArray3& cell) const;
  //void findGradeIndex(const long64& cell_key,
  //    IntArray3& cell) const;


  //   void applyVelocityBC(BodySP& body) const;

  //   bool intersection(const Point3D& point, const Vector3D& ray,
  //                     Point3D& hitPoint) const;

private:

  IntVector3D  d_num_particles_per_cell;   // Bryan's ppe
  Point3D      d_lower;                    //         X0
  Point3D      d_upper;                    //         X1
  IntVector3D  d_node_counts;              //         Nc
  Vector3D     d_num_ghost;                //         nGhost
  Vector3D     d_cell_size;                //         dX
  double       d_time;                     //         t
  double       d_time_final;               //         tf
  double       d_delT;                     //         dt
  int          d_iteration;                //         it
  double       d_tol;                      //         tol

  /* TODO: Hooman implement MPMBC (called in MPM2D) */
  // MPMBCs d_vel_BC;              //         bcs
  MPMShapeFunctionP d_shape;

};  // end class

}  // end namespace
#endif
