#ifndef EMU2DC_BOND_H
#define EMU2DC_BOND_H

#include <NodeP.h>

namespace Emu2DC {

  class Bond 
  {
  public: 
    Bond();
    Bond(NodeP node1, NodeP node2);
    virtual ~Bond();

    void first(const NodeP node);
    void second(const NodeP node);
    void isBroken(const bool broken);

    NodeP first() const;
    NodeP second() const;
    bool isBroken() const;
    
    bool operator==(const Bond& bond) const;

  private:
    NodeP d_node1;
    NodeP d_node2;
    bool d_broken;
  };  // end class

} // end namespace

#endif
