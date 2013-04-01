#ifndef MATITI_BONDMATERIAL_H
#define MATITI_BONDMATERIAL_H

#include <vector>

namespace Matiti {

  class BondMaterial
  {
    public:
      BondMaterial();
      ~BondMaterial();

    private:

      int d_matType;
      double d_density;
      double d_youngModulus;
      double d_strainEnergy;
      double d_damageIndex;

      // prevent blank creation and copying
      BondMaterial(const BondMaterial& bond);

  };
} // End namespace Matiti

#endif // MATITI_BONDMATERIAL_H
