#include <CCA/Components/MPM/ParticleCreator/ParticleGeometryHelpers.h>
#include <Core/GeometryPiece/BoxGeometryPiece.h>
#include <Core/GeometryPiece/CylinderGeometryPiece.h>
#include <Core/GeometryPiece/SphereGeometryPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Patch.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <cmath>

namespace Vaango::Helpers {

void createCylinderSurfacePoints(const Uintah::CylinderGeometryPiece* cyl,
                                 const Uintah::Patch* patch,
                                 std::vector<Uintah::Point>& points) {
  Uintah::Point bottom = cyl->bottom();
  Uintah::Point top    = cyl->top();
  double rad   = cyl->radius();
  Vector axis  = top - bottom;
  double h     = axis.length();
  axis.normalize();

  // Find two vectors orthogonal to the axis
  Uintah::Vector v1, v2;
  axis.findOrthogonal(v1, v2);

  const int n_theta = 16;
  const int n_z     = 8;
  
  // Generate points along the cylindrical surface
  for (int i_z = 0; i_z <= n_z; ++i_z) {
    double z = static_cast<double>(i_z) / static_cast<double>(n_z) * h;
    Point p_z = bottom + axis * z;
    
    for (int i_t = 0; i_t < n_theta; ++i_t) {
      double theta = static_cast<double>(i_t) / static_cast<double>(n_theta) * 2.0 * M_PI;
      Point p = p_z + (v1 * cos(theta) + v2 * sin(theta)) * rad;
      
      if (patch->containsPoint(p)) {
        points.push_back(p);
      }
    }
  }
  
  // Generate points on the caps (top and bottom)
  const int n_r = 3;
  for (int i_r = 1; i_r < n_r; ++i_r) {
    double r = static_cast<double>(i_r) / static_cast<double>(n_r) * rad;
    
    for (int i_t = 0; i_t < n_theta; ++i_t) {
      double theta = static_cast<double>(i_t) / static_cast<double>(n_theta) * 2.0 * M_PI;
      Vector radial = (v1 * cos(theta) + v2 * sin(theta)) * r;
      
      Point p_bot = bottom + radial;
      Point p_top = top + radial;
      
      if (patch->containsPoint(p_bot)) {
        points.push_back(p_bot);
      }
      if (patch->containsPoint(p_top)) {
        points.push_back(p_top);
      }
    }
  }
}

void createSphereSurfacePoints(const Uintah::SphereGeometryPiece* sphere,
                               const Uintah::Patch* patch,
                               std::vector<Uintah::Point>& points) {
  Uintah::Point center = sphere->origin();
  double rad   = sphere->radius();
  
  const int n_theta = 16;
  const int n_phi   = 8;
  
  // Generate points using spherical coordinates (excluding poles)
  for (int i_p = 1; i_p < n_phi; ++i_p) {
    double phi = static_cast<double>(i_p) / static_cast<double>(n_phi) * M_PI;
    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    
    for (int i_t = 0; i_t < n_theta; ++i_t) {
      double theta = static_cast<double>(i_t) / static_cast<double>(n_theta) * 2.0 * M_PI;
      
      Uintah::Point p = center + Uintah::Vector(rad * sin_phi * cos(theta),
                                rad * sin_phi * sin(theta),
                                rad * cos_phi);
      
      if (patch->containsPoint(p)) {
        points.push_back(p);
      }
    }
  }
  
  // Add the poles
  Uintah::Point p_north = center + Uintah::Vector(0, 0, rad);
  Uintah::Point p_south = center + Uintah::Vector(0, 0, -rad);
  
  if (patch->containsPoint(p_north)) {
    points.push_back(p_north);
  }
  if (patch->containsPoint(p_south)) {
    points.push_back(p_south);
  }
}

void
createBoxSurfacePoints(const Uintah::BoxGeometryPiece* box_piece,
                       const Uintah::Patch* patch,
                       std::vector<Uintah::Point>& points)
{
  Uintah::Box box    = box_piece->getBoundingBox();
  Uintah::Point low  = box.lower();
  Uintah::Point high = box.upper();
  auto cell_spacing = patch->dCell();

  const int n_x = (high.x() - low.x()) / cell_spacing.x();
  const int n_y = (high.y() - low.y()) / cell_spacing.y();
  const int n_z = (high.z() - low.z()) / cell_spacing.z();

  auto add_point = [&](const Point& p) {
    if (patch->containsPoint(p)) {
      points.push_back(p);
    }
  };

  // Faces: x = low.x and x = high.x
  for (int i = 0; i <= n_y; ++i) {
    double y = low.y() + (high.y() - low.y()) * i / n_y;
    for (int j = 0; j <= n_z; ++j) {
      double z = low.z() + (high.z() - low.z()) * j / n_z;
      add_point(Point(low.x(), y, z));
      add_point(Point(high.x(), y, z));
    }
  }

  // Faces: y = low.y and y = high.y
  for (int i = 0; i <= n_x; ++i) {
    double x = low.x() + (high.x() - low.x()) * i / n_x;
    for (int j = 0; j <= n_z; ++j) {
      double z = low.z() + (high.z() - low.z()) * j / n_z;
      add_point(Point(x, low.y(), z));
      add_point(Point(x, high.y(), z));
    }
  }

  // Faces: z = low.z and z = high.z
  for (int i = 0; i <= n_x; ++i) {
    double x = low.x() + (high.x() - low.x()) * i / n_x;
    for (int j = 0; j <= n_y; ++j) {
      double y = low.y() + (high.y() - low.y()) * j / n_y;
      add_point(Point(x, y, low.z()));
      add_point(Point(x, y, high.z()));
    }
  }
}

void createBoundingBoxCornerPoints(const Uintah::Box& box,
                                   const Uintah::Point& center,
                                   const Uintah::Patch* patch,
                                   std::vector<Uintah::Point>& points) {
  // Define all 8 corners of the bounding box
  Uintah::Point corners[8];
  corners[0] = box.lower();
  corners[1] = Uintah::Point(box.upper().x(), box.lower().y(), box.lower().z());
  corners[2] = Uintah::Point(box.upper().x(), box.upper().y(), box.lower().z());
  corners[3] = Uintah::Point(box.lower().x(), box.upper().y(), box.lower().z());
  corners[4] = Uintah::Point(box.lower().x(), box.lower().y(), box.upper().z());
  corners[5] = Uintah::Point(box.upper().x(), box.lower().y(), box.upper().z());
  corners[6] = box.upper();
  corners[7] = Uintah::Point(box.lower().x(), box.upper().y(), box.upper().z());

  Uintah::Vector eps(1e-10, 1e-10, 1e-10);
  
  for (int k = 0; k < 8; ++k) {
    Point p = corners[k];
    
    // Adjust corner point slightly toward center to avoid boundary issues
    Vector to_center = center - p;
    if (to_center.length() > 1e-9) {
      to_center.normalize();
      p = p + to_center * 1e-9;
    }

    // Check if point or slightly adjusted point is in patch
    if (patch->containsPoint(p) || patch->containsPoint(p - eps)) {
      // Ensure uniqueness (important for 2D where corners may overlap)
      bool exists = false;
      for (const auto& existing_p : points) {
        if ((existing_p - p).length2() < 1.0e-16) {
          exists = true;
          break;
        }
      }
      
      if (!exists) {
        points.push_back(p);
      }
    }
  }
}

Uintah::Point adjustPointToPatch(const Uintah::Point& point, const Uintah::Patch* patch) {
  if (patch->containsPoint(point)) {
    return point;
  }
  
  Uintah::Vector eps(1e-10, 1e-10, 1e-10);
  
  // Try small adjustments to get the point inside the patch
  if (patch->containsPoint(point - eps)) {
    return point - eps * 0.1;
  }
  
  if (patch->containsPoint(point + eps)) {
    return point + eps * 0.1;
  }
  
  return point;
}

} // namespace Helpers
