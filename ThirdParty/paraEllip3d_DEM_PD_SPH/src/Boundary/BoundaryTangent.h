#ifndef BOUNDARY_TANGENT_H
#define BOUNDARY_TANGENT_H

#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>

namespace dem {

/////////////////////////////////////
class BoundaryTangent
{
public:
  std::size_t particleId;
  Vec tangentForce;
  Vec tangentDisp;
  bool tangentLoading;
  Vec tangentDispStart;
  REAL tangentPeak;

public:
  BoundaryTangent()
    : particleId(0)
    , tangentForce(0)
    , tangentDisp(0)
    , tangentLoading(false)
    , tangentDispStart(0)
    , tangentPeak(0)
  {
  }

  BoundaryTangent(std::size_t _particleId, Vec _v1, Vec _v2, bool _b, Vec _v3,
                  REAL _tp)
    : particleId(_particleId)
    , tangentForce(std::move(_v1))
    , tangentDisp(std::move(_v2))
    , tangentLoading(_b)
    , tangentDispStart(std::move(_v3))
    , tangentPeak(_tp)
  {
  }
};

} // end namespace dem

#endif
