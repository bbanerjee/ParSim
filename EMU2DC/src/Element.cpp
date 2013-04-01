#include <Element.h>

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
