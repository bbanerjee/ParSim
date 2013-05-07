#ifndef EMU2DC_ELEMENT_H
#define EMU2DC_ELEMENT_H

#include <NodeP.h>
#include <NodePArray.h>
#include <Geometry/Point3D.h>

namespace Emu2DC {
    
  class Element {
      
  public:

    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::Element& elem);

  public:

    Element();
    Element(const int& id, const NodePArray& nodes);
    ~Element();

    void initialize(const int id, const NodePArray& nodes);

    inline int id() const {return d_id;}
    const NodePArray& nodes() const {return d_nodes;}
    inline int numNodes() const {return d_nodes.size();}

    void computeGeometry2D(double& area, double& xlength, double& ylength) const;

    void computeVolume();
    void computeVolume2D();
    void computeVolume3D();

    double volume() const {return d_volume;}

  protected:

    int d_id;
    int d_dimension;
    double d_volume;
    NodePArray d_nodes; 

  private:

    double computeVolumeTetrahedron(const Point3D& p0, const Point3D& p1, const Point3D& p2,
                                    const Point3D& p3) const;

    // Prevent copy construction and operator=
    Element(const Element& element);
    Element& operator=(const Element& element);

  }; // end class

} // end namespace

#endif


