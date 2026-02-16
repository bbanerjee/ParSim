#ifndef PARTICLE_GEOMETRY_HELPERS_H
#define PARTICLE_GEOMETRY_HELPERS_H

#include <vector>

// Forward declare to avoid loops
namespace Uintah {
class CylinderGeometryPiece;
class SphereGeometryPiece;
class BoxGeometryPiece;
class Patch;
class Point;
class Box;
} // namespace Uintah

namespace Vaango::Helpers {

/**
 * Creates surface visualization points for a cylinder geometry.
 * Generates points on the cylindrical surface, top cap, and bottom cap.
 *
 * @param cyl The cylinder geometry piece
 * @param patch The patch to check for point containment
 * @param points Output vector to append generated points
 */
void
createCylinderSurfacePoints(const Uintah::CylinderGeometryPiece* cyl,
                            const Uintah::Patch* patch,
                            std::vector<Uintah::Point>& points);

/**
 * Creates surface visualization points for a sphere geometry.
 * Generates points distributed across the sphere surface using spherical
 * coordinates.
 *
 * @param sphere The sphere geometry piece
 * @param patch The patch to check for point containment
 * @param points Output vector to append generated points
 */
void
createSphereSurfacePoints(const Uintah::SphereGeometryPiece* sphere,
                          const Uintah::Patch* patch,
                          std::vector<Uintah::Point>& points);

/**
 * Creates surface visualization points for a box geometry.
 * Generates points on all six faces of the box.
 *
 * @param box_piece The box geometry piece
 * @param patch The patch to check for point containment
 * @param points Output vector to append generated points
 */
void
createBoxSurfacePoints(const Uintah::BoxGeometryPiece* box_piece,
                       const Uintah::Patch* patch,
                       std::vector<Uintah::Point>& points);

/**
 * Creates corner points for a bounding box, slightly adjusted toward the
 * center. This ensures the object is visible to neighbor patches.
 *
 * @param box The bounding box
 * @param center The center point to adjust corners toward
 * @param patch The patch to check for point containment
 * @param points Output vector to append generated corner points
 */
void
createBoundingBoxCornerPoints(const Uintah::Box& box,
                              const Uintah::Point& center,
                              const Uintah::Patch* patch,
                              std::vector<Uintah::Point>& points);

/**
 * Attempts to adjust a point to be inside the patch if it's very close to the
 * boundary.
 *
 * @param point The point to potentially adjust
 * @param patch The patch to check containment
 * @return The adjusted point if successful, original point otherwise
 */
Uintah::Point
adjustPointToPatch(const Uintah::Point& point, const Uintah::Patch* patch);

} // namespace Vaango::Helpers

#endif // PARTICLE_GEOMETRY_HELPERS_H
