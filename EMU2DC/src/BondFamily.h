#ifndef EMU2DC_BONDFAMILY_H
#define EMU2DC_BONDFAMILY_H

//****************************************************************************
//*  Purpose    : This module includes subroutines to calculate the      *
//*        bond-family structures                    *
//****************************************************************************

namespace Emu2DC {

  class BondFamily {  

  public:

    BondFamily();
    ~BondFamily();

    void getFamilyReference(const Node* node,
                            NodeArray& family);

    void getFamilyDeformed(const Node* node,
                           NodeArray& family);

    void sortNodesReference();
    void sortNodesDeformed();

  private:

    // prevent copying
    BondFamily(const BondFamily& family);
    BondFamily& operator=(const BondFamily& family);

  };  // end class
}  // End namespace 
#endif 

