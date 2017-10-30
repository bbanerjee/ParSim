#include <Core/Const/Constants.h>
#include <Core/Geometry/Box.h>
#include <Core/Geometry/OrientedBox.h>
#include <iostream>

namespace dem {

OrientedBox::OrientedBox(const Box& box) {
  d_center = box.getCenter();
  d_axes[0] = Vec(1, 0, 0);
  d_axes[1] = Vec(0, 1, 0);
  d_axes[2] = Vec(0, 0, 1);
  d_half_len[0] = box.getDimx()*0.5;
  d_half_len[1] = box.getDimy()*0.5;
  d_half_len[2] = box.getDimz()*0.5;
}

void 
OrientedBox::normalize_axes() {
  for (auto& axis: d_axes) {
    axis.normalizeInPlace();
  }
}

std::vector<Vec> 
OrientedBox::vertices() const
{
  std::vector<Vec> vertices;  vertices.reserve(8);
  std::array<Vec, 3> axes = {{d_axes[0]*d_half_len[0],
                              d_axes[1]*d_half_len[1],
                              d_axes[2]*d_half_len[2]}};
  vertices.push_back(d_center - axes[0] - axes[1] - axes[2]);
  vertices.push_back(d_center + axes[0] - axes[1] - axes[2]);
  vertices.push_back(d_center + axes[0] + axes[1] - axes[2]);
  vertices.push_back(d_center - axes[0] + axes[1] - axes[2]);
  vertices.push_back(d_center - axes[0] - axes[1] + axes[2]);
  vertices.push_back(d_center + axes[0] - axes[1] + axes[2]);
  vertices.push_back(d_center + axes[0] + axes[1] + axes[2]);
  vertices.push_back(d_center - axes[0] + axes[1] + axes[2]);
  return vertices;                  
}

bool 
OrientedBox::containsPoint(const Vec& point) const
{
  Vec ww = d_center - point;
  REAL dist_a = dot(ww, d_axes[0]);
  REAL dist_b = dot(ww, d_axes[1]);
  REAL dist_c = dot(ww, d_axes[2]);
  if (dist_a > d_half_len[0] || dist_b > d_half_len[1] || dist_c > d_half_len[2]) {
    return false;
  }
  return true;
}

/**
 * Algorithm from : Eberly, 2008, 
 *                  "Dynamic Collision Detection using Oriented Bounding Boxes"
 * https://www.geometrictools.com/
 */
bool 
OrientedBox::intersects(const OrientedBox& box) const
{
  // Create the D vector
  Vec D = box.d_center - d_center;

  // Arrays to store A_i . D and B_j . D as they are computed
  std::array<REAL, 3> AD, BD;

  // Storage for the C matrix: C_ij = A_i . B_j as it is computed
  Matrix3 C(0.0);

  // Compute R0, R1, R2 for L = A0, A1, A2
  for (auto ii = 0u; ii < 3; ++ii) {
    REAL R0 = d_half_len[ii];

    REAL R1 = 0.0;
    for (auto jj = 0u; jj < 3; ++jj) {
      C(ii, jj) = dot(d_axes[ii], box.d_axes[jj]);
      R1 += box.d_half_len[jj] * std::abs(C(ii,jj));
    }

    AD[ii] = dot(d_axes[ii], D);
    REAL R = std::abs(AD[ii]);

    if (R > R0 + R1) return false;
  }

  // Compute R0, R1, R2 for L = B0, B1, B2
  for (auto jj = 0u; jj < 3; ++jj) {

    REAL R0 = 0.0;
    for (auto ii = 0u; ii < 3; ++ii) {
      R0 += d_half_len[ii] * std::abs(C(ii,jj));
    }

    REAL R1 = box.d_half_len[jj];

    BD[jj] = dot(box.d_axes[jj], D);
    REAL R = std::abs(BD[jj]);

    if (R > R0 + R1) return false;
  }

  // Compute R0, R1, R2 for L = A_i x (B0, B1, B2)
  std::array<std::array<int, 2>, 3> indices = {{{{1, 2}}, {{2, 0}}, {{0, 1}}}};
  for (auto ii = 0u; ii < 3; ++ii) {
    auto i1 = indices[ii][0];
    auto i2 = indices[ii][1];
    for (auto jj = 0u; jj < 3; ++jj) {
      auto j1 = indices[jj][0];
      auto j2 = indices[jj][1];

      REAL R0 = d_half_len[i1] * std::abs(C(i2, jj)) +
                d_half_len[i2] * std::abs(C(i1, jj));

      REAL R1 = box.d_half_len[j1] * std::abs(C(ii, j2)) +
                box.d_half_len[j2] * std::abs(C(ii, j1));

      REAL R = std::abs(C(i1, jj) * AD[i2] - C(i2, jj) * AD[i1]);

      if (R > R0 + R1) return false;
    }
  }

  return true; 
}

void 
OrientedBox::rotate(REAL angle, const Vec& axis)
{
  Matrix3 rot(angle, axis);
  d_axes[0] = rot*d_axes[0];
  d_axes[1] = rot*d_axes[1];
  d_axes[2] = rot*d_axes[2];
}

void 
OrientedBox::translate(const Vec& dist)
{
  d_center += dist;
}

std::ostream&
operator<<(std::ostream& os, const OrientedBox& b)
{
  os << " center = " << b.d_center 
     << " axis_a = " << b.d_axes[0]
     << " axis_b = " << b.d_axes[1]
     << " axis_c = " << b.d_axes[2]
     << " half_len_a = " << b.d_half_len[0]
     << " half_len_b = " << b.d_half_len[1]
     << " half_len_c = " << b.d_half_len[2]
     << "\n";
  return os;
}

} // namespace dem
