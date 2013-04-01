#ifndef EMU2DC_DOMAIN_H
#define EMU2DC_DOMAIN_H

#include <array>

namespace Emu2DC {

  class Domain {

  public:  
    typedef std::array<double, 3> Array3;
    typedef std::array<int, 3> IntArray3;

  public:  

    Domain() ;
    ~Domain();

    Domain(const Array3& lower, const Array3& upper);

    Domain(const Array3& lower, const Array3& upper, const IntArray3& numcells);
    
    Domain(const Array3& lower, const Array3& upper, const double& horizon);

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

    IntArray3 d_numcells;

    // Don't allow copy
    Domain(const Domain& dom);
    Domain& operator=(const Domain& dom);

  };  // end class
}  // end namespace
#endif
