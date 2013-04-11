#include <Element.h>
#include <Node.h>
#include <Types.h>
#include <algorithm>
#include <cmath>

using namespace Emu2DC;

// The element constructor
Element::Element()
 : d_id(0), d_nodes(0)
{
}
  
// The element destructor
Element::~Element()
{
}

void 
Element::initialize(const int id, const NodePArray& nodes)
{
  d_id = id;
  d_nodes = nodes;
}

void Element::computeGeometry2D(double& area, double& xlength, double& ylength) const
{
  // Set up coordinate array
  std::vector<double> xcoords, ycoords;
  for (constNodePIterator iter = d_nodes.begin(); iter != d_nodes.end(); iter++) {
    NodeP node = *iter;
    Array3 pos = node->position();
    xcoords.push_back(pos[0]);
    ycoords.push_back(pos[1]);
  }
  xcoords.push_back(xcoords[0]);
  ycoords.push_back(ycoords[0]);

  // Compute area
  area = 0.0;
  for (unsigned int ii = 0; ii < xcoords.size()-1; ii++) {
    area += (xcoords[ii]*ycoords[ii+1]-xcoords[ii+1]*ycoords[ii]);
  }
  area = std::abs(0.5*area);

  // Compute xlength, ylength
  double xmax = *(std::max_element(xcoords.begin(), xcoords.end()));
  double xmin = *(std::min_element(xcoords.begin(), xcoords.end()));
  double ymax = *(std::max_element(ycoords.begin(), ycoords.end()));
  double ymin = *(std::min_element(ycoords.begin(), ycoords.end()));
  xlength = std::abs(xmax-xmin);
  ylength = std::abs(ymax-ymin);
}
