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
  Vec tgtForce;
  Vec tgtDisp;
  bool tgtLoading;
  Vec tgtDispStart;
  REAL tgtPeak;

public:
  BoundaryTangent()
    : particleId(0)
    , tgtForce(0)
    , tgtDisp(0)
    , tgtLoading(false)
    , tgtDispStart(0)
    , tgtPeak(0)
  {
  }

  BoundaryTangent(std::size_t _particleId, Vec _v1, Vec _v2, bool _b, Vec _v3,
                  REAL _tp)
    : particleId(_particleId)
    , tgtForce(std::move(_v1))
    , tgtDisp(std::move(_v2))
    , tgtLoading(_b)
    , tgtDispStart(std::move(_v3))
    , tgtPeak(_tp)
  {
  }
};

} // end namespace dem

#endif
