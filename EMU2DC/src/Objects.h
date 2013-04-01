#ifndef EMU2DC_OBJECTS
#define EMU2DC_OBJECTS

namespace Emu2DC {
    
  // This structure defines the node type
  class Node {

    public:

      int id;
      int mat_type;
      double horizon_size;
      bool iflag;  // iflag = 1: node is a hanging node; 
                   // iflag = 0: node is not a hanging node

      double volume;
      double density;
      double young;

      double strain_energy = 0;
      double damage_index;

      double* pos;  // array 
      double* disp;  // array
      double* veloc;  // array
      double* accel;  // array
      double* new_veloc;  // array
      double* new_disp;  // array
      double* old_disp  // array
      double* force;  // array

      int* bc;
    
      int nnodeelements;
      int nodeelements[10]; // The id of the elements adjacent to this node, 
                            // can be 1,2,3 or 4 elements
      double nodeelements_size[10]; // The size(area) of the elements adjacent to this node, 
                                    // can be 1,2,3 or 4 elements
      int nodeelements_depth[10]; // The depth of each element adjacent to this node
      int nodeelements_nhanging_nodes[10]; // The number of hanging nodes for each element 
                                           // connected to this node
  };

  //   Node position on an element
  //    _______
  //  1|       |2 
  //   |       |
  //   |       |
  //  4|_______|3
  //      
  class Element {
      
    public:

      Element();
      ~Element();

      // This function is created to assign one element to another. 
      void Element_attribution(const Element* element);

      
    public:

      int id;

      Node* elementnodes; // Array
    
      Node* node1;
      Node* node2;
      Node* node3;
      Node* node4;

      int depth;
      bool leath = true;    // tell us if this element is a leaf of the quadtree stucture
      bool root = false;    
      bool dif_level_refine = false; // flag that tell us if the element should be refined because 
                                     // a dif. level of refinment >2 situation happened
      bool strain_energy_refine;

      Element* child1;     // we can improve this by changing this pointer to point to a id number 
                           // instead of a data structure  
      Element* child2;     // this id number is from the global array of elements
      Element* child3;
      Element* child4;

      Element* father;
      Element* quad_element;

      int n_neighbors;
      int* neighborhood; // Array

    private:

      // Prevent copy construction and operator=
      Element(const Element* element);
      Element& operator=(const Element* element);

  };

      
  class Line {
    
    public:

      double x1, y1, x2, y2;
    
  };

// The element constructor
Element::Element()
{
  d_id = 0;
  d_elementnodes = 0;
  d_n_neighbors = 0;
  d_neighborbood = 0;
  d_node1 = 0;
  d_node2 = 0;
  d_node3 = 0;
  d_node4 = 0;
  d_depth = 0;
  d_leath = true;
  d_root = false;
  d_dif_level_refine = false;

  d_child1 = 0;
  d_child2 = 0;
  d_child3 = 0;
  d_child4 = 0;
  d_father = 0;
  d_quad_element = 0;
}
  
// The element destructor
Element::~Element()
{
  delete d_elementnodes;
  delete d_neighborhood;
}

// This function is created to assign one element to another. 
void
Element::Element_attribution(const Element* element)
{
  d_id = element->d_id;

  d_elementnodes = new Node[4];
  for (int ii = 0; ii < 4; ii++) {
    d_elementnodes[ii] = element->d_elementnodes[ii];
  }

  d_n_neighbors = element->d_n_neighbors;
  d_neighborhood = new int[d_n_neighbors];
  for (int ii =0; ii < d_n_neighbors; ii++) {
    d_neighborhood[ii] = element->d_neighborhood[ii];
  }

  d_node1 = element->d_node1;
  d_node2 = element->d_node2; 
  d_node3 = element->d_node3; 
  d_node4 = element->d_node4; 

  d_depth = element->d_depth;
  d_leath = element->d_leath;
  d_root = element->d_root;
  d_child1 = element->d_child1;
  d_child2 = element->d_child2;
  d_child3 = element->d_child3;
  d_child4 = element->d_child4;
  d_father = element->d_father;

  d_quad_element => element->d_quad_element;
}
