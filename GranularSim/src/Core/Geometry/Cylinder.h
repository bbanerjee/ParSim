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
    : d_radius(0)
    , d_height(0)
    , d_center(0)
  {
  }

  Cylinder(REAL r, REAL h, Vec c)
    : d_radius(r)
    , d_height(h)
    , d_center(std::move(c))
  {
  }

  Cylinder(const Cylinder& cy)
  {
    d_radius = cy.d_radius;
    d_height = cy.d_height;
    d_center = cy.d_center;
  }

  REAL radius() const { return d_radius; }
  REAL getHeight() const { return d_height; }
  REAL volume() const { return Pi * d_radius * d_radius * d_height; }
  Vec center() const { return d_center; }

  void setRadius(REAL r) { d_radius = r; }
  void setHeight(REAL h) { d_height = h; }
  void setCenter(Vec v) { d_center = v; }

  Vec randomPoint() const;
  void print() const;

private:
  REAL d_radius;
  REAL d_height;
  Vec d_center;
};
}

#endif
