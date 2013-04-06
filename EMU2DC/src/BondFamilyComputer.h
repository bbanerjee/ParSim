#ifndef EMU2DC_BONDFAMILYCOMPUTER_H
#define EMU2DC_BONDFAMILYCOMPUTER_H

#include <Domain.h>
#include <Node.h>
#include <NodeP.h>

//************************************** 
/** 
  * @file
  * @author  Biswajit Banerjee (2013)
  * @version 1.0
  *
  * @section LICENSE
  *
  * Copyright (c) 2013 Callaghan Innovation
  *
  * @section DESCRIPTION
  *
  * This class provides methods for computing the bond-family structures
  * for a node.
  */

namespace Emu2DC {

  class BondFamilyComputer {  

  public:

    /**
     * Create a BondfamilyComputer object and the associated Cell-NodeP map
     *
     * @param domain Reference to the domain object
     * @param nodeList Reference to the vector of NodeP objects inside the domain
     */
    BondFamilyComputer(const Domain& domain,
                       const NodePArray& nodeList);

    ~BondFamilyComputer();

    void getInitialFamily(const NodeP& node,
                          NodePArray& family);

    void getCurrentFamily(const NodeP& node,
                          NodePArray& family);

    void sortNodesReference();
    void sortNodesDeformed();

  private:

    // Prevent empty construction
    BondFamilyComputer();

    // prevent copying
    BondFamilyComputer(const BondFamilyComputer& family);
    BondFamilyComputer& operator=(const BondFamilyComputer& family);

  };  // end class
}  // End namespace 
#endif 

