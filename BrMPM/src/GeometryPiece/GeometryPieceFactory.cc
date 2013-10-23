#include <MPMDatawarehouse.h>
#include <GeometryPiece/GeometryPieceFactory.h>
#include <GeometryPiece/GeometryPiece.h>
#include <GeometryPiece/BoxGeometryPiece.h>
//#include <GeometryPiece/GeometryReader.h>
#include <Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

using namespace BrMPM;

GeometryPiece* 
GeometryPieceFactory::create(Uintah::ProblemSpecP& ps,
                             MPMDatawarehouseP& dw)
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
    //return (new BoxGeometryPiece(geom_ps, dw));
    return NULL;
  } else if (geom_type == "file") {
    //return (new GeometryReader(geom_ps, dw));
    return NULL;
  } else {
    std::ostringstream out;
    out << "**ERROR** Unknown geometry type" << geom_type;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }
}
