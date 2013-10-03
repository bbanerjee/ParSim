#ifndef __MATITI_FAMILYCOMPUTER_H__
#define __MATITI_FAMILYCOMPUTER_H__

#include <Domain.h>
#include <Node.h>
#include <NodeP.h>
#include <NodePArray.h>
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
  * This class provides methods for computing the node-family structures
  * for a node.
  */

namespace Matiti {

  class FamilyComputer {  

  public:

    /**
     * Create an empty BondfamilyComputer object
     */
    FamilyComputer();

    ~FamilyComputer();

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
                          const Domain& domain,
                          NodePArray& family) const;

    /**
     *  Finds the family of a node: all the nodes inside the horizon of the node 
     *    The family is based on the current nodal positions
     *
     * @param node shared_ptr to the node object
     * @param family Reference to the vector of NodeP objects that makes up the family of node
     */
    void getCurrentFamily(NodeP node, 
                          const Domain& domain,
                          NodePArray& family) const;

  private:

    // Store the cell-node map
    CellNodePMap d_map;

    // prevent copying
    FamilyComputer(const FamilyComputer& family);
    FamilyComputer& operator=(const FamilyComputer& family);

  };  // end class
}  // End namespace 
#endif 

