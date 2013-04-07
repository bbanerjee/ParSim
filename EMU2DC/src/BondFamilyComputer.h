#ifndef EMU2DC_BONDFAMILYCOMPUTER_H
#define EMU2DC_BONDFAMILYCOMPUTER_H

#include <Domain.h>
#include <Node.h>
#include <NodeP.h>
#include <CellNodePMap.h>

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
     * Create an empty BondfamilyComputer object
     */
    BondFamilyComputer();

    ~BondFamilyComputer();

    /**
     *  Find which cells the nodes sit in and create a unordered map that maps nodes to cells
     *
     * @param domain Reference to the domain object
     * @param nodeList Reference to the vector of NodeP objects inside the domain
     */
    void createCellNodeMap(const Domain& domain,
                           const NodePArray& nodeList);

    /**
     *  Clear existing map and create a new map
     *
     * @param domain Reference to the domain object
     * @param nodeList Reference to the vector of NodeP objects inside the domain
     */
    void updateCellNodeMap(const Domain& domain,
                           const NodePArray& nodeList);

    /**
     *  Print cell-node map
     */
    void printCellNodeMap() const;
    void printCellNodeMap(const IntArray3& cell) const;

    /**
     *  Finds the family of a node: all the nodes inside the horizon of the node 
     *    The family is based on the initial nodal positions
     *
     * @param node shared_ptr to the node object
     * @param family Reference to the vector of NodeP objects that makes up the family of node
     */
    void getInitialFamily(NodeP node,
                          NodePArray& family) const;

    void getCurrentFamily(NodeP node,
                          NodePArray& family) const;

    void sortNodesReference();
    void sortNodesDeformed();

  private:

    // Store the cell-node map
    CellNodePMap d_map;

    // prevent copying
    BondFamilyComputer(const BondFamilyComputer& family);
    BondFamilyComputer& operator=(const BondFamilyComputer& family);

  };  // end class
}  // End namespace 
#endif 

