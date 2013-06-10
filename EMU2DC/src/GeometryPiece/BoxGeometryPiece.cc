#include <GeometryPiece/BoxGeometryPiece.h>
#include <Exception.h>
#include <Node.h>
#include <Element.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Geometry/Vector.h>

#include <iostream>

using namespace Emu2DC; 

BoxGeometryPiece::BoxGeometryPiece(Uintah::ProblemSpecP& ps,
                                   NodePArray& nodes,
                                   ElementPArray& elements)
{
  d_name = "box";
  Uintah::Vector lower, upper;
  ps->require("lower", lower);
  ps->require("upper", upper);
  d_box = Box3D(Point3D(lower[0],lower[1],lower[2]), Point3D(upper[0],upper[1],upper[2]));
  if (d_box.isDegenerate()) {
    std::ostringstream out;
    out << "**ERROR** The box geometry piece is degenerate. Lower = " << d_box.lower()
        << " Upper = " << d_box.upper() << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }
  ps->require("num_edge_nodes", d_num_edge_nodes);
  if (d_num_edge_nodes < 1) d_num_edge_nodes = 1;
}

BoxGeometryPiece::~BoxGeometryPiece()
{
}

Box3D 
BoxGeometryPiece::boundingBox() const
{
  return d_box;
}


bool 
BoxGeometryPiece::inside (const Point3D& pt) const
{
  return d_box.contains(pt); 
}


std::string 
BoxGeometryPiece::name() const
{
  return d_name;
}

void
BoxGeometryPiece::createNodes(NodePArray& nodes, ElementPArray& elements)
{
}

