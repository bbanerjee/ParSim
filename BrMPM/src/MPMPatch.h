#ifndef __MATITI_MPMPATCH_H__
#define __MATITI_MPMPATCH_H__

#include <Domain.h>
#include <Types.h>
#include <BodySP.h>
#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>
#include <Geometry/IntVector3D.h>
<<<<<<< HEAD
#include <VelocityBCSPArray.h>
=======

>>>>>>> 508431166e437d17f4913ab5b0c99723ffbc80f6
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <iostream>

namespace BrMPM {

  class MPMPatch
 {

  public:  

//    friend std::ostream& operator<<(std::ostream& out, const Matiti::Domain& domain);

  public:  

    MPMPatch() ;
     ~MPMPatch();

    MPMPatch(const Point3D& lower, const Point3D& upper);
    MPMPatch(const Point3D& lower, const Point3D& upper, const IntArray3& numCells);
    MPMPatch(const Uintah::ProblemSpecP& ps);

 //   Domain(const Point3D& lower, const Point3D& upper, const IntArray3& numCells);
    
 //   Domain(const Point3D& lower, const Point3D& upper, const double& horizon);

    void initialize(const Uintah::ProblemSpecP& ps);

<<<<<<< HEAD
  
    const int& ghost() const {return d_num_ghost;}
=======
    Vector3D& nGhost() {return d_num_ghost;}
    IntVector3D& nC()  {return d_node_counts;}
    Vector3D& dX()  {return d_cell_size;}
    Point3D& x0() {return d_lower;}

    const Vector3D& ghost() const {return d_num_ghost;}
>>>>>>> 508431166e437d17f4913ab5b0c99723ffbc80f6
 //   const double& thick() const {return d_thick;}
    const int& particlesperelement() const {return d_num_particles_per_cell;}
   // const IntArray3& numGrids() const;
    const double totalGrids() const;

    const Point3D lower() const {return d_lower;}
    const Point3D upper() const {return d_upper;}


    const Vector3D& cellSize()  {return d_cell_size;}
    const Vector3D& numGrids()  {return d_num_grids;}
    const std::vector<Point3D> gridsPosition()  {return d_gridsPosition;}

    void findGradeIndex(const Point3D& point,
                       IntArray3& cell) const;
    void findGradeIndex(const long64& cell_key,
                       IntArray3& cell) const;

    bool insidePatch(const Point3D& point) const;
    bool allInsidePatch(const std::vector<Point3D> points) const;

 //   void applyVelocityBC(BodySP& body) const;

 //   bool intersection(const Point3D& point, const Vector3D& ray,
 //                     Point3D& hitPoint) const;

  private:
<<<<<<< HEAD
   
    Point3D d_lower;    // X0
    Point3D d_upper;    // X1
    Vector3D d_node_counts;  // Nc
    Vector3D d_cellsize;     // dX
    IntVector3D d_num_cells;
    
   // int d_t_initial;
   // int d_t_final;
    int d_num_ghost;        //nGhost
=======

    Vector3D d_num_ghost;      // Bryan's nGhost
    Point3D d_lower;           //         X0
    IntVector3D d_node_counts; //         nC
    Vector3D d_cell_size;       //         dX

    //IntVector3D d_nC;
    //Vector3D d_dX;
    //Point3D d_X0;

 //   Point3D d_lower;
 //   Point3D d_upper;
   
    Point3D d_upper;
    
   // int d_t_initial;
   // int d_t_final;
>>>>>>> 508431166e437d17f4913ab5b0c99723ffbc80f6
    double d_tol;                    //tolerance
   // double d_thick;
    int d_num_particles_per_cell;

    
    double xcoord, ycoord, zcoord;
    

    Vector3D d_num_grids;
    Point3D d_grids;

    std::vector<Point3D>  d_gridsPosition;  
 //   VelocityBCSPArray d_vel_BC;
   
  };  // end class
}  // end namespace
#endif
