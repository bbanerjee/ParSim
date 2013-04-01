#ifndef EMU2DC_BOND_H
#define EMU2DC_BOND_H

#include <Node.h>

namespace Emu2DC {

  class Bond 
  {
  public: 
    Bond();
    Bond(Node* node1, Node* node2);
    virtual ~Bond();

    void first(const Node* node);
    void second(const Node* node);
    void isBroken(const bool broken);

    Node* first() const;
    Node* second() const;
    bool isBroken() const;
    
    bool operator==(const Bond& bond) const;

  private:
    Node* d_node1;
    Node* d_node2;
    bool d_broken;
  };  // end class

} // end namespace

#endif
