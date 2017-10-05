#ifndef CYLINDER_H
#define CYLINDER_H

#include <Core/Const/Constants.h>
#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>

namespace dem {

class Cylinder
{

public:
  Cylinder()
    : radius(0)
    , height(0)
    , center(0)
  {
  }

  Cylinder(REAL r, REAL h, Vec c)
    : radius(r)
    , height(h)
    , center(std::move(c))
  {
  }

  Cylinder(const Cylinder& cy)
  {
    radius = cy.radius;
    height = cy.height;
    center = cy.center;
  }

  REAL getRadius() const { return radius; }
  REAL getHeight() const { return height; }
  REAL getVolume() const { return Pi * radius * radius * height; }
  Vec getCenter() const { return center; }

  void setRadius(REAL r) { radius = r; }
  void setHeight(REAL h) { height = h; }
  void setCenter(Vec v) { center = v; }

  Vec randomPoint() const;
  void print() const;

private:
  REAL radius;
  REAL height;
  Vec center;
};
}

#endif
