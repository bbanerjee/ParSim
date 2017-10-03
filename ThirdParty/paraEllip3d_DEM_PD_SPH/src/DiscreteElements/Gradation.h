#ifndef GRADATION_H
#define GRADATION_H

#include <Core/Types/realtypes.h>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <InputOutput/IOUtils.h>
#include <InputOutput/zenxml/xml.h>
#include <cstddef>
#include <vector>
#include <iostream>
#include <fstream>

namespace dem {

class Gradation
{

public:
  Gradation()
    : sieveNum()
    , percent()
    , size()
    , ratioBA()
    , ratioCA()
  {
  }

  Gradation(std::size_t sn, std::vector<REAL> v1, std::vector<REAL> v2, REAL ba,
            REAL ca)
    : sieveNum(sn)
    , percent(v1)
    , size(v2)
    , ratioBA(ba)
    , ratioCA(ca)
  {
  }

  inline void set(std::size_t sn, std::vector<REAL> v1, std::vector<REAL> v2,
                  REAL ba, REAL ca)
  {
    sieveNum = sn;
    percent = v1;
    size = v2;
    ratioBA = ba;
    ratioCA = ca;
  }

  void initializeFromCSVFile(std::ifstream& ifs)
  {
    std::size_t sieveNum;
    ifs >> sieveNum;
    std::vector<REAL> percent(sieveNum), size(sieveNum);
    for (auto i = 0u; i < sieveNum; ++i)
      ifs >> percent[i] >> size[i];
    REAL ratio_ba, ratio_ca;
    ifs >> ratio_ba >> ratio_ca;
    set(sieveNum, percent, size, ratio_ba, ratio_ca);
  }

  // **TODO** Add validity checks
  void initializeFromXMLFile(zen::XmlIn& ps)
  {
    auto sieve_ps = ps["Sieves"];
    if (sieve_ps) {
      std::size_t numSieves;
      sieve_ps.attribute("number", numSieves);

      std::string percentPassingStr;
      sieve_ps["percent_passing"](percentPassingStr);
      std::vector<REAL> percentPassing = Ellip3D::Util::convertStrArray<REAL>(percentPassingStr);

      std::string sizeStr;
      sieve_ps["size"](sizeStr);
      std::vector<REAL> size = Ellip3D::Util::convertStrArray<REAL>(sizeStr);

      REAL ratio_ba, ratio_ca;
      sieve_ps["sieve_ratio"]["ratio_ba"](ratio_ba);
      sieve_ps["sieve_ratio"]["ratio_ca"](ratio_ca);

      set(numSieves, percentPassing, size, ratio_ba, ratio_ca);
    }
  }

  std::size_t getSieveNum() const { return sieveNum; }
  std::vector<REAL>& getPercent() { return percent; }
  const std::vector<REAL>& getPercent() const { return percent; }
  std::vector<REAL>& getSize() { return size; }
  const std::vector<REAL>& getSize() const { return size; }
  REAL getPtclRatioBA() const { return ratioBA; }
  REAL getPtclRatioCA() const { return ratioCA; }
  void setPtclRatioBA(REAL ba) { ratioBA = ba; }
  void setPtclRatioCA(REAL ca) { ratioCA = ca; }

  REAL getPtclMaxRadius() const { return size[0]; }
  REAL getPtclMinRadius() const { return size[sieveNum - 1] * ratioCA; }

private:
  std::size_t sieveNum; // sieveNum == percent.size() == size.size()
  std::vector<REAL> percent;
  std::vector<REAL> size;
  REAL ratioBA; // ratio of radius b to radius a
  REAL ratioCA; // ratio of radius c to radius a

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& sieveNum;
    ar& percent;
    ar& size;
    ar& ratioBA;
    ar& ratioCA;
  }
};

} // namespace dem

#endif
