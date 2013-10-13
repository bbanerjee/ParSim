#ifndef __MATITI_MPMPATCH_H__
#define __MATITI_MPMPATCH_H__

#include <Domain.h>
#include <Types.h>
#include <BodySP.h>
#include <Geometry/Point3D.h>
#include <VelocityBCSPArray.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <iostream>

namespace BrMPM {

  class MPMPatch : public Domain
 {

  public:  

//    friend std::ostream& operator<<(std::ostream& out, const Matiti::Domain& domain);

  public:  

    MPMPatch() ;
    virtual ~MPMPatch();

    MPMPatch(const Uintah::ProblemSpecP& ps);

 //   Domain(const Point3D& lower, const Point3D& upper, const IntArray3& numCells);
    
 //   Domain(const Point3D& lower, const Point3D& upper, const double& horizon);

    void initialize(const Uintah::ProblemSpecP& ps);

  
    const int& ghost() const {return d_ghost;}
    const double& thick() const {return d_thick;}
    const int& particlesperelement() const {return d_particlesperelement;}
    const IntArray3& numGrids() const;
    const double totalGrids() const;
    const std::vector<Point3D> gridsPosition() const {return d_gridsPosition;}

    void findGradeIndex(const Point3D& point,
                       IntArray3& cell) const;
    void findGradeIndex(const long64& cell_key,
                       IntArray3& cell) const;

    bool insidePatch(const Point3D& point) const;

 //   void applyVelocityBC(BodySP& body) const;

 //   bool intersection(const Point3D& point, const Vector3D& ray,
 //                     Point3D& hitPoint) const;

  private:

 //   Point3D d_lower;
 //   Point3D d_upper;
   
    Point3D d_lower;
    Point3D d_upper;
    Vector3D d_node_counts;
    
    int d_t_initial;
    int d_t_final;
    int d_ghost;
    double d_tol                    //tolerance
    double d_thick;
    int d_particlesperelement;

    double xcoord, ycoord, zcoord;
    

    IntArray3 d_num_grids;
    Point3D d_grids;

    std::vector<Point3D>  d_gridsPosition;  
 //   VelocityBCSPArray d_vel_BC;
   
  };  // end class
}  // end namespace
#endif
