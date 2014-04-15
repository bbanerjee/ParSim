#ifndef __MATITI_DENSITY_H__
#define __MATITI_DENSITY_H__

#include <Pointers/NodeP.h>
#include <Pointers/DensitySP.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <vector>

namespace Matiti {

  
  class Density {

  public:

    Density();
    Density(const Density& den);
    virtual ~Density();

    void clone(const DensitySP& den);

    void initialize(Uintah::ProblemSpecP& ps);

    void nodeDensity (const NodeP& node, double& node_density);

    inline const double& ringWidth() const { return d_ring_width; }
    inline void ringWidth(const double& width) { d_ring_width = width; }

  protected:

    double remind (double lengthPeriod, double nodePos);

    double density (const std::vector<double>& polyCoeff, double reminder);

  private:

    double d_ring_width;
    std::vector<double> d_poly_coeffs;

 }; // end class

}; // end namespace

#endif

