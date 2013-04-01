#ifndef MATITI_BONDDAMAGEMODEL_H
#define MATITI_BONDDAMAGEMODEL_H

namespace Matiti {

  class BondDamageModel
  {
    public:
      BondDamageModel();
      ~BondDamageModel();

      void initializeBrokenBonds();
      void computeCriticalStrain();
      void updateBrokenBonds();

    private:

      // prevent copying
      BondDamageModel(const BondDamageModel& dm);
      BondDamageModel& operator=(const BondDamageModel& dm);

  };
} // End namespace Matiti

#endif // MATITI_BONDDAMAGEMODEL_H
