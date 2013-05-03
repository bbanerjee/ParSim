#ifndef EMU2DC_ELEMENT_H
#define EMU2DC_ELEMENT_H

#include <NodeP.h>
#include <NodePArray.h>

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

    void computeGeometry2D(double& area, double& xlength, double& ylength) const;

  protected:

    int d_id;
    NodePArray d_nodes; 

  private:

    // Prevent copy construction and operator=
    Element(const Element& element);
    Element& operator=(const Element& element);

  }; // end class

} // end namespace

#endif


