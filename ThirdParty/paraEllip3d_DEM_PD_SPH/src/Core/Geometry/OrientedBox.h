#ifndef __ELLIP3D_CORE_GEOM_ORIENTED_BOX_H__
#define __ELLIP3D_CORE_GEOM_ORIENTED_BOX_H__

#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>
#include <vector>

namespace dem {

class OrientedBox
{
private:

  Vec d_center;
  Vec d_axis_a;
  Vec d_axis_b;
  Vec d_axis_c;
  REAL d_half_len_a;
  REAL d_half_len_b;
  REAL d_half_len_c;

  void normalize_axes();

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_center;
    ar& d_axis_a;
    ar& d_axis_b;
    ar& d_axis_c;
    ar& d_half_len_a;
    ar& d_half_len_b;
    ar& d_half_len_c;
  }

public:

  OrientedBox(const Vec& center, 
              const Vec& ax_a, const Vec& ax_b, const Vec& ax_c,
              REAL rad_a, REAL rad_b, REAL rad_c)
    : d_center(center)
    , d_axis_a(ax_a)
    , d_axis_b(ax_b)
    , d_axis_c(ax_c)
    , d_half_len_a(rad_a)
    , d_half_len_b(rad_b)
    , d_half_len_c(rad_c)
  {
    normalize_axes();
  }

  std::vector<Vec> vertices() const;

  bool containsPoint(const Vec& pt) const;

  friend std::ostream& operator<<(std::ostream& os, const OrientedBox& box);

};

} // namespace dem

#endif
