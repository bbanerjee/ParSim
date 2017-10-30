#ifndef __ELLIP3D_CORE_GEOM_ORIENTED_BOX_H__
#define __ELLIP3D_CORE_GEOM_ORIENTED_BOX_H__

#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>

#include <array>
#include <vector>

namespace dem {

class Box;

class OrientedBox
{
private:

  Vec d_center;
  std::array<Vec, 3> d_axes; // axis_a, axis_b, axis_c
  std::array<REAL, 3> d_half_len; // len_a, len_b, len_c

  void normalize_axes();

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_center;
    ar& d_axes;
    ar& d_half_len;
  }

public:

  OrientedBox(const Vec& center, 
              const Vec& ax_a, const Vec& ax_b, const Vec& ax_c,
              REAL rad_a, REAL rad_b, REAL rad_c)
    : d_center(center)
    , d_axes({{ax_a, ax_b, ax_c}})
    , d_half_len({{rad_a, rad_b, rad_c}})
  {
    normalize_axes();
  }

  OrientedBox(const Box& box);

  const Vec& center() const {return d_center;}
  const Vec& axis(int idx) const {return d_axes.at(idx);}
  REAL extent(int idx) const {return d_half_len.at(idx);}

  std::vector<Vec> vertices() const;

  bool containsPoint(const Vec& pt) const;

  bool intersects(const OrientedBox& box) const;

  void rotate(REAL angle, const Vec& axis);
  void translate(const Vec& dist);

  friend std::ostream& operator<<(std::ostream& os, const OrientedBox& box);

};

} // namespace dem

#endif
