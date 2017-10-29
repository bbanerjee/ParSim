#include <Core/Const/Constants.h>
#include <Core/Geometry/OrientedBox.h>
#include <iostream>

namespace dem {

void 
OrientedBox::normalize_axes() {
  d_axis_a.normalizeInPlace();
  d_axis_b.normalizeInPlace();
  d_axis_c.normalizeInPlace();
}

std::vector<Vec> 
OrientedBox::vertices() const
{
  std::vector<Vec> vertices;  vertices.reserve(8);
  Vec aa = d_axis_a*d_half_len_a;
  Vec bb = d_axis_b*d_half_len_b;
  Vec cc = d_axis_c*d_half_len_c;
  vertices.push_back(d_center - aa - bb - cc);
  vertices.push_back(d_center + aa - bb - cc);
  vertices.push_back(d_center + aa + bb - cc);
  vertices.push_back(d_center - aa + bb - cc);
  vertices.push_back(d_center - aa - bb + cc);
  vertices.push_back(d_center + aa - bb + cc);
  vertices.push_back(d_center + aa + bb + cc);
  vertices.push_back(d_center - aa + bb + cc);
  return vertices;                  
}

bool 
OrientedBox::containsPoint(const Vec& point) const
{
  Vec ww = d_center - point;
  REAL dist_a = dot(ww, d_axis_a);
  REAL dist_b = dot(ww, d_axis_b);
  REAL dist_c = dot(ww, d_axis_c);
  if (dist_a > d_half_len_a || dist_b > d_half_len_b || dist_c > d_half_len_c) {
    return false;
  }
  return true;
}

std::ostream&
operator<<(std::ostream& os, const OrientedBox& b)
{
  os << " center = " << b.d_center 
     << " axis_a = " << b.d_axis_a 
     << " axis_b = " << b.d_axis_b 
     << " axis_c = " << b.d_axis_c 
     << " half_len_a = " << b.d_half_len_a 
     << " half_len_b = " << b.d_half_len_b 
     << " half_len_c = " << b.d_half_len_c 
     << "\n";
  return os;
}

} // namespace dem
