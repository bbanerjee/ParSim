#ifndef GRADATION_H
#define GRADATION_H

#include <Core/Types/realtypes.h>
#include <vector>
#include <cstddef>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

namespace dem {

class Gradation {

public:
  Gradation() : sieveNum(), percent(), size(), ratioBA(), ratioCA() {}

  Gradation(std::size_t sn, std::vector<REAL> v1, std::vector<REAL> v2, REAL ba,
            REAL ca)
      : sieveNum(sn), percent(v1), size(v2), ratioBA(ba), ratioCA(ca) {}

  std::size_t getSieveNum() const { return sieveNum; }
  std::vector<REAL> &getPercent() { return percent; }
  const std::vector<REAL> &getPercent() const { return percent; }
  std::vector<REAL> &getSize() { return size; }
  const std::vector<REAL> &getSize() const { return size; }
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
  void serialize(Archive &ar, const unsigned int version) {
    ar &sieveNum;
    ar &percent;
    ar &size;
    ar &ratioBA;
    ar &ratioCA;
  }

};

} // namespace dem

#endif
