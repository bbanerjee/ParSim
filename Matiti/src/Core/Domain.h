#ifndef __MATITI_DOMAIN_H__
#define __MATITI_DOMAIN_H__

#include <Types/Types.h>
#include <Pointers/BodySP.h>
#include <Geometry/Point3D.h>
#include <Containers/VelocityBCSPArray.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Core/Geometry/Vector.h>
#include <iostream>

namespace Matiti {

  class Domain {

  public:  

    friend std::ostream& operator<<(std::ostream& out, const Matiti::Domain& domain);

  public:  

    Domain() ;
    virtual ~Domain();

    Domain(const Point3D& lower, const Point3D& upper);

    Domain(const Point3D& lower, const Point3D& upper, const IntArray3& numCells);
    
    Domain(const Point3D& lower, const Point3D& upper, const SCIRun::Vector& cellSize);

    virtual void initialize(const Uintah::ProblemSpecP& ps);

    const Point3D& lower() const;
    const Point3D& upper() const;
    const double& xrange() const;
    const double& yrange() const;
    const double& zrange() const;
    const IntArray3& numCells() const;
    const SCIRun::Vector& cellSize() const;
    const double totalCells() const;

    void findCellIndex(const Point3D& point,
                       IntArray3& cell) const;
    void findCellIndex(const long64& cell_key,
                       IntArray3& cell) const;

    bool inside(const Point3D& point) const;

    void applyVelocityBC(BodySP& body) const;

    bool intersection(const Point3D& point, const Vector3D& ray,
                      Point3D& hitPoint) const;

  private:

    Point3D d_lower;
    Point3D d_upper;

    double d_xrange;
    double d_yrange;
    double d_zrange;

    IntArray3 d_num_cells;
    SCIRun::Vector d_cell_size;

    VelocityBCSPArray d_vel_BC;

    // Don't allow copy
    Domain(const Domain& dom);
    Domain& operator=(const Domain& dom);

  };  // end class
}  // end namespace
#endif
