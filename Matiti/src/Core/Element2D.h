#ifndef MATITI_ELEMENT_2D_H
#define MATITI_ELEMENT_2D_H

#include <Pointers/NodeP.h>
#include <Containers/NodePArray.h>
#include <Geometry/Point3D.h>

/*
 Special class for computing surface areas
 */
namespace Matiti {
    
  class Element2D {
      
  public:

    friend std::ostream& operator<<(std::ostream& out, const Matiti::Element2D& elem);

  public:

    Element2D();
    Element2D(const NodeP& node1, const NodeP& node2, const NodeP& node3);
    Element2D(const NodeP& node1, const NodeP& node2, const NodeP& node3, 
              const NodeP& node4);
    Element2D(const NodePArray& nodes);
    ~Element2D();

    void initialize(const NodePArray& nodes);

    const NodePArray& nodes() const {return d_nodes;}
    inline int numNodes() const {return d_nodes.size();}

    bool hasNode(const NodeP& node) const;
    bool isSubset(const NodePArray& nodes) const;

    void computeArea();

    double area() const {return d_area;}
    double surfaceArea() const {return d_surface_area;}

  private:

    double computeAreaTriangle(const Point3D& p1, const Point3D& p2, 
                               const Point3D& p3) const;

    double computeAreaQuadrilateral(const Point3D& p1, const Point3D& p2, 
                                    const Point3D& p3, const Point3D& p4) const;
  private:

    double d_area;
    double d_surface_area; // area = 0 if not on surface
    NodePArray d_nodes; 

    // Prevent copy construction and operator=
    Element2D(const Element2D& element);
    Element2D& operator=(const Element2D& element);

  }; // end class

} // end namespace

#endif


