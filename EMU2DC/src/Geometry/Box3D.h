#ifndef __EMU2DC_BOX3D_H__
#define __EMU2DC_BOX3D_H__

#include <Geometry/Point3D.h>

#include <iostream>

namespace Emu2DC {

  class Box3D 
  {
  public:

    friend std::ostream& operator<<(std::ostream& os, const Box3D& box);

  public:

    Box3D();
    Box3D(const Box3D& box);
    Box3D(const Point3D& lower, const Point3D& upper);

    ~Box3D();

    Box3D& operator=(const Box3D& box);

    bool overlaps(const Box3D& box, double epsilon=1.0e-6) const;

    bool contains(const Point3D& pt) const;

    bool isDegenerate();

    Point3D lower() const;
    Point3D upper() const;

  private:

    void correctBoundingBox();

  private:

    Point3D d_lower;
    Point3D d_upper;

  }; // end class Box3D
 
} // End namespace 

#endif //ifndef __EMU2DC_BOX3D_H__
