#ifndef __EMU2DC_DOMAIN_H__
#define __EMU2DC_DOMAIN_H__

#include <Types.h>
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

    Domain(const Array3& lower, const Array3& upper);

    Domain(const Array3& lower, const Array3& upper, const IntArray3& numCells);
    
    Domain(const Array3& lower, const Array3& upper, const double& horizon);

    void initialize(const Uintah::ProblemSpecP& ps);

    const Array3& lower() const;
    const Array3& upper() const;
    const double& horizon() const;
    const double& xrange() const;
    const double& yrange() const;
    const double& zrange() const;
    const IntArray3& numCells() const;
    const double totalCells() const;

    void findCellIndex(const Array3& point,
                       IntArray3& cell) const;

    bool inside(const Array3& point) const;

  private:

    Array3 d_lower;
    Array3 d_upper;

    double d_horizon;
    double d_xrange;
    double d_yrange;
    double d_zrange;

    IntArray3 d_num_cells;
    VelocityBCSPArray d_vel_BC;

    // Don't allow copy
    Domain(const Domain& dom);
    Domain& operator=(const Domain& dom);

  };  // end class
}  // end namespace
#endif
