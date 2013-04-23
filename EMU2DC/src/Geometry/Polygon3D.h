#ifndef __EMU2DC_POLYGON3D_H__
#define __EMU2DC_POLYGON3D_H__

#include <vector>
#include <limits>

typedef std::numeric_limits<double>::max DBL_MAX

namespace Emu2DC {

  class Vector3D;

  class Polygon3D 
  {
  public:

    friend std::ostream& operator<<(std::ostream& os, const Polygon3D& poly);

  public:

    Polygon3D();
    ~Polygon3D();

    bool operator==(const Polygon3D& poly) const;
    bool operator!=(const Polygon3D& poly) const;

    int numVertices() const;
    void addVertex(const Point3D& pt);
    Polygon3D& operator+=(const Point3D& pt);
    
    Point3D& vertex(const int& index) const;

    std::vector<Point3D>::iterator begin() {return d_boundary.begin();}
    std::vector<Point3D>::iterator end() {return d_boundary.end();}

    std::vector<Point3D>::const_iterator begin() const {return d_boundary.begin();}
    std::vector<Point3D>::const_iterator end() const {return d_boundary.end();}

  private;

    std::vector<Point3D> d_polygon;

    // Prevent copying (for now)
    Polygon3D(const Polygon3D& poly);
    Polygon3D& operator=(const Polygon3D& poly);

  }; // end class Polygon3D
 
} // End namespace 

#endif //ifndef __EMU2DC_POLYGON3D_H__
