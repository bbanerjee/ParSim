#ifndef EMU2DC_DATAWAREHOUSE_H
#define EMU2DC_DATAWAREHOUSE_H

#include <Node.h>
#include <Element.h>
#include <Line.h>
#include <string>
#include <vector>
#include <array>

namespace Emu2DC {

  class DataWarehouse {

  public:

    typedef std::vector<Node*> NodeArray;
    typedef std::vector<Node*>::iterator NodeArrayIter;
    typedef std::vector<Element*> ElemArray;
    typedef std::vector<Element*>::iterator ElemArrayIter;

  public:
    
    DataWarehouse();
    ~DataWarehouse();

    inline void getNodeList(NodeArray& nodes) const
    {
      nodes = d_nodes;
    }

    inline void setNodeList(const NodeArray& nodes) 
    {
      d_nodes = nodes;
    }

    inline void getOriginalNodeList(NodeArray& nodes) const
    {
      nodes = d_original_nodes;
    }

    inline void setOriginalNodeList(const NodeArray& nodes) 
    {
      d_original_nodes = nodes;
    }

    inline void getElementList(ElemArray& elems) const
    {
      elems = d_elements;
    }

    inline void setElementList(const ElemArray& elems) 
    {
      d_elements = elems;
    }

    inline void getOriginalElementList(ElemArray& elems) const
    {
      elems = d_old_elements;
    }

    inline void setOriginalElementList(const ElemArray& elems) 
    {
      d_old_elements = elems;
    }

  private:

     NodeArray d_nodes;            // global vector of nodes
     NodeArray d_original_nodes;   

     ElemArray d_elements;                   // global vector of  leath elements
     ElemArray d_old_elements;               // old global vector of leath elements
     ElemArray d_quadtree_elements;          // initial vector of elements( input elements)
     ElemArray d_original_quadtree_elements; 

     //_____________________ Solver Global variables _____________________
     std::vector<std::vector<int> > bond;
     std::vector<double> critical_strain;

     std::vector<bool> mine;

     std::vector<std::vector<int> > broke;
     std::vector<bool > omitt;

     std::vector<double> edt, edt_add;
     std::vector<double> ut1new, ut2new, ut1, ut2;
     std::vector<double> ut1old, ut2old;
     std::vector<double> f1, f2;
     std::vector<double> wt;

     std::vector<double> acc1, acc2;  // Acceleration of the nodes: acc1 = acceleration in the x-direction, 
                                      // acc2 = acceleration in the y-direction 

     std::vector<double> lc1bd, lc2bd, u1bd, u2bd, v1bd, v2bd, tendbd;

     std::vector<double> vis;
     std::vector<int> node_type, nofail;

     std::vector<double> fcoefm;       // coefficents of bond
     std::vector<int> fnorm;           // normalization of bond force
     std::vector<double> crit_exten;   // critical extension of bonds
     std::vector<double> damage_index; // Each node has a value for damage_index, which is the ratio of broken bonds 
                                       // for a node to the total number of bonds for the same node in the reference 
                                       // configuration
     std::vector<double> volnod;       // volume of node
     std::vector<int> nodbd, mbdtyp;
     std::vector<double> spsum;

     std::vector<int> family;
     std::vector<int> nodes_in_bin;
     std::vector<int> def_family;
     std::vector<int> nodes_in_bin_def;

     std::vector<double> radnod;

   
} // end namespace
#endif
