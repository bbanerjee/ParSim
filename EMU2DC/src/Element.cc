#include <Element.h>
#include <algorithm>

// The element constructor
Element::Element()
{
  d_id = 0;
  d_elementnodes = 0;
  d_n_neighbors = 0;
  d_neighborbood = 0;
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

void Element::computeGeometry(double& area, double& xlength, double& ylength) const
{
  // Set up coordinate array
  Array3 pos = {0.0, 0.0, 0.0}; 
  std::vector<double> xcoords, ycoords;
  std::vector<Node*>::iterator iter = d_elementnodes.begin();
  for (; iter != d_elementnodes.end(); ++iter) {
    d_elementnodes[*iter]->getPosition(pos);
    xcoords.push_back(pos[0]);
    ycoords.push_back(pos[1]);
  }
  iter = d_elementnodes.begin();
  d_elementnodes[*iter]->getPosition(pos);
  xcoords.push_back(pos[0]);
  ycoords.push_back(pos[1]);

  // Compute area
  area = 0.0;
  for (int ii = 0; ii < xcoords.size()-1; ii++) {
    area += (xcoords[ii]*ycoords[ii+1]-xcoords[ii+1]*ycoords[ii]);
  }
  area = std::fabs(0.5*area);

  // Compute xlength, ylength
  double xmax = *(std::max_element(xcoords.begin(), xcoords.end()));
  double xmin = *(std::min_element(xcoords.begin(), xcoords.end()));
  double ymax = *(std::max_element(ycoords.begin(), ycoords.end()));
  double ymin = *(std::min_element(ycoords.begin(), ycoords.end()));
  xlength = std:fabs(xmax-xmin);
  ylength = std:fabs(ymax-ymin);

}
