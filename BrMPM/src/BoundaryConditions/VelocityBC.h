#ifndef __MATITI_VELOCITYBC_H__
#define __MATITI_VELOCITYBC_H__

#include <NodeP.h>
#include <GeometryMath/Point3D.h>
#include <GeometryMath/Vector3D.h>
#include <GeometryMath/Polygon3D.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <vector>

namespace BrMPM {
  
  class VelocityBC {
  
  public:

    enum class BCType {
      Symmetry=0,
      Periodic=1,
      Wall=2,
      Outlet=3
    };

    enum class FaceType {
      Xminus=0,
      Xplus=1,
      Yminus=2,
      Yplus=3,
      Zminus=4,
      Zplus=5
    };

    VelocityBC();
    virtual ~VelocityBC();

    void initialize(Uintah::ProblemSpecP& ps);
    void apply(NodeP& node, 
               const Point3D& hitPoint,
               const Point3D& domain_min, 
               const Point3D& domain_max) const;

    const BCType& bcType() const {return d_bc;}
    const FaceType& face() const {return d_face;}
    inline double getRestitution() const {return d_restitution;}

  private:

    void updateVelocityAndPosition(NodeP& node, 
                                   const Point3D& hitPoint,
                                   const Vector3D& normal) const; 

  private:

    // Coefficient of restitution
    double d_restitution;

    // Applied velocity BC
    BCType d_bc;

    // Face on which BC is applied
    FaceType d_face;

    // Points that describe the regions of application of velocity BC on the domain
    // boundary. 
    std::vector<BCType> d_area_bc;
    std::vector<Polygon3D> d_area_boundary;

    // prevent copying
    VelocityBC(const VelocityBC& bc);
    VelocityBC& operator=(const VelocityBC& bc);

  }; // end class

} // end namespace
#endif

