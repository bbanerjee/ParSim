#ifndef EMU2DC_BODY_H
#define EMU2DC_BODY_H

#include <Material.h>
#include <NodePArray.h>
#include <ElementPArray.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>

namespace Emu2DC {

  class Body
  {
  public:
   
    Body();
    ~Body();

    void initialize(ProblemSpecP& ps);

    // **WARNING** One mat for now.  A body can have more than one material.
    inline int matID() const {return d_mat_id;}
    const NodePArray& nodes() const {return d_nodes;}
    const ElementPArray& elements() const {return d_elements;}

  protected:

    void readNodeFile(const std::string& fileName);

  private:

    int d_mat_id;
    NodePArray d_nodes;
    ElementPArray d_elements;

  };
} // end namespace

#endif
