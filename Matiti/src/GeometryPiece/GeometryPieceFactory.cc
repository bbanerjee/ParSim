#include <GeometryPiece/GeometryPieceFactory.h>
#include <GeometryPiece/GeometryPiece.h>
#include <GeometryPiece/BoxGeometryPiece.h>
#include <GeometryPiece/GeometryReader.h>
#include <Core/Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

using namespace Matiti;

GeometryPiece* 
GeometryPieceFactory::create(Uintah::ProblemSpecP& ps,
                             NodePArray& nodes,
                             ElementPArray& elements)
{
  // Get the geometry 
  Uintah::ProblemSpecP geom_ps = ps->findBlock("Geometry");
  if (!geom_ps) {
    std::ostringstream out;
    out << "**ERROR** No geometry information found for body";
    throw Exception(out.str(), __FILE__, __LINE__); 
  }
  
  // Get the attribute of the geometry
  std::string geom_type;
  if (!geom_ps->getAttribute("type", geom_type)) {
    std::ostringstream out;
    out << "**ERROR** Geometry does not have type information" << geom_type;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }

  // Create geometry
  if (geom_type == "box") {
    return (new BoxGeometryPiece(geom_ps, nodes, elements));
  } else if (geom_type == "file") {
    return (new GeometryReader(geom_ps, nodes, elements));
  } else {
    std::ostringstream out;
    out << "**ERROR** Unknown geometry type" << geom_type;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }
}
