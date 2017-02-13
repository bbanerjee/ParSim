#include <Core/Math/IntVec.h>
#include <Core/Geometry/Box.h>
#include <Boundary/Containers.h>
#include <Boundary/BoundaryReader.h>
#include <Boundary/PlaneBoundary.h>
#include <Boundary/CylinderBoundary.h>
#include <InputOutput/zenxml/xml.h>

using namespace dem;

void
BoundaryReader::read(const char* str,
                     Box& container,
                     Box& grid,
                     BoundaryPArray& boundaries) const
{
  std::ifstream ifs(str);
  if (!ifs) {
    std::cout << "**ERROR**: Could not read boundary information from "
              << str << " in " << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  REAL x1, y1, z1, x2, y2, z2;
  ifs >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
  container.set(x1, y1, z1, x2, y2, z2);
  // compute grid assumed to be the same as container, change in
  // scatterParticle() if necessary.
  grid.set(x1, y1, z1, x2, y2, z2);

  boundaries.clear();
  std::size_t boundaryNum;
  std::size_t type;
  ifs >> boundaryNum;
  for (std::size_t i = 0; i < boundaryNum; ++i) {
    ifs >> type;
    if (type == 1) // plane boundary
      boundaries.push_back(std::make_shared<PlaneBoundary>(type, ifs));
    else if (type == 2) // cylindrical boundary
      boundaries.push_back(std::make_shared<CylinderBoundary>(type, ifs));
  }

  ifs.close();
}

bool
BoundaryReader::readXML(const std::string& inputFileName,
                        Box& container,
                        Box& grid,
                        BoundaryPArray& boundaries) const
{
  // Read the input file
  zen::XmlDoc doc;
  try {
    std::cout << "Input file name= " << inputFileName << "\n";
    doc = zen::load(inputFileName);
  } catch (const zen::XmlFileError& err) {
    std::cout << "*ERROR** Could not read input file " << inputFileName << "\n";
    std::cout << "    Error # = " << err.lastError << "\n";
    return false;
  } catch (const zen::XmlParsingError& err) {
    std::cout << "*ERROR** Could not read input file " << inputFileName << "\n";
    std::cout << "    Parse Error in line: " << err.row + 1
              << " col: " << err.col << "\n";
    return false;
  }

  // Check whether this is the right type of input file
  if (doc.root().getNameAs<std::string>() != "Ellip3D_input") {
    std::cout << "*ERROR** Could not find tag <Ellip3D_input> in input file "
              << inputFileName << "\n";
    return false;
  }

  // Load the document into input proxy for easier element access
  zen::XmlIn ps(doc);

  // Read the title
  std::string title;
  if (!ps["Meta"]["title"](title)) {
    std::cout << "*ERROR** Could not find boundary title in input file "
              << inputFileName << "\n";
    std::cout << "  Add the <title> tag inside a <Meta> tag\n";
    return false;
  }
  std::cout << "title = " << title << "\n";

  // Read the boundary information
  auto boundary_ps = ps["Boundary"];
  if (!boundary_ps) {
    std::cout << "**ERROR** <Boundary> tag not found. \n";
    return false;
  }

  // Read the container dimensions
  std::string vecStr;
  if (!boundary_ps["containerMin"](vecStr)) {
    std::cout << "**ERROR** Container min. position not found in boundary geometry\n";
    std::cout << "  Add the <containerMin> [x, y, z] </containerMin> tag.";
    return false;
  }
  Vec boxMin = Vec::fromString(vecStr);

  if (!boundary_ps["containerMax"](vecStr)) {
    std::cout << "**ERROR** Container max. position not found in boundary geometry\n";
    std::cout << "  Add the <containerMax> [x, y, z] </containerMax> tag.";
    return false;
  }
  Vec boxMax = Vec::fromString(vecStr);

  container.set(boxMin.getX(), boxMin.getY(), boxMin.getZ(),
                boxMax.getX(), boxMax.getY(), boxMax.getZ());

  // compute grid assumed to be the same as container, change in
  // scatterParticle() if necessary.
  grid.set(boxMin.getX(), boxMin.getY(), boxMin.getZ(),
           boxMax.getX(), boxMax.getY(), boxMax.getZ());

  BoundaryType type; 
  boundaries.clear();
  for (auto bound_ps = ps["boundary"]; bound_ps; bound_ps.next()) {
    std::string boundaryType;
    bound_ps.attribute("type", boundaryType);
    BoundaryId id;
    bound_ps.attribute("id", id);
    switch (getEnum(boundaryType)) {
      case PLANE:
        type = 1; 
        boundaries.push_back(std::make_shared<PlaneBoundary>(type, id, bound_ps));
        break;
      case CYLINDER:
        type = 2; 
        boundaries.push_back(std::make_shared<CylinderBoundary>(type, id, bound_ps));
        break;
      case NONE:
        break;
    }
  }

  return true;
}
