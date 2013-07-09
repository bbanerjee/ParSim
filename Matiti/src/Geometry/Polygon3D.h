#ifndef __MATITI_POLYGON3D_H__
#define __MATITI_POLYGON3D_H__

#include <Geometry/Point3D.h>

#include <vector>

namespace Matiti {

  class Polygon3D 
  {
  public:

    friend std::ostream& operator<<(std::ostream& os, const Polygon3D& poly);

  public:

    Polygon3D();
    Polygon3D(const std::vector<Point3D>& poly);
    Polygon3D(const Polygon3D& poly);
    ~Polygon3D();

    bool operator==(const Polygon3D& poly) const;

    unsigned int numVertices() const;
    void addVertex(const Point3D& pt);
    Polygon3D& operator+=(const Point3D& pt);
    
    const Point3D& vertex(const int& index) const;
    const Point3D& operator[](const int& index) const;

    std::vector<Point3D>::iterator begin() {return d_vertices.begin();}
    std::vector<Point3D>::iterator end() {return d_vertices.end();}

    std::vector<Point3D>::const_iterator begin() const {return d_vertices.begin();}
    std::vector<Point3D>::const_iterator end() const {return d_vertices.end();}

  private:

    std::vector<Point3D> d_vertices;

    // Prevent copying (for now)
    Polygon3D& operator=(const Polygon3D& poly);

  }; // end class Polygon3D
 
} // End namespace 

#endif //ifndef __MATITI_POLYGON3D_H__
