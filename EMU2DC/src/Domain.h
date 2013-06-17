#ifndef __EMU2DC_DOMAIN_H__
#define __EMU2DC_DOMAIN_H__

#include <Types.h>
#include <Geometry/Point3D.h>
#include <VelocityBCSPArray.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <iostream>

namespace Emu2DC {

  class Domain {

  public:  

    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::Domain& domain);

  public:  

    Domain() ;
    ~Domain();

    Domain(const Point3D& lower, const Point3D& upper);

    Domain(const Point3D& lower, const Point3D& upper, const IntArray3& numCells);
    
    Domain(const Point3D& lower, const Point3D& upper, const double& horizon);

    void initialize(const Uintah::ProblemSpecP& ps);

    const Point3D& lower() const;
    const Point3D& upper() const;
    const double& horizon() const;
    const double& xrange() const;
    const double& yrange() const;
    const double& zrange() const;
    const IntArray3& numCells() const;
    const double totalCells() const;

    void findCellIndex(const Point3D& point,
                       IntArray3& cell) const;
    void findCellIndex(const long64& cell_key,
                       IntArray3& cell) const;

    bool inside(const Point3D& point) const;

  private:

    Point3D d_lower;
    Point3D d_upper;

    double d_xrange;
    double d_yrange;
    double d_zrange;
    double d_horizon;

    IntArray3 d_num_cells;
    VelocityBCSPArray d_vel_BC;

    // Don't allow copy
    Domain(const Domain& dom);
    Domain& operator=(const Domain& dom);

  };  // end class
}  // end namespace
#endif
