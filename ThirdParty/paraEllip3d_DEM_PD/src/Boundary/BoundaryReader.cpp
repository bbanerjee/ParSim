#include <Boundary/BoundaryReader.h>
#include <Boundary/Containers.h>
#include <Boundary/CylinderBoundary.h>
#include <Boundary/PlaneBoundary.h>
#include <Core/Geometry/Box.h>
#include <Core/Math/IntVec.h>
#include <InputOutput/json/json.hpp>
#include <InputOutput/zenxml/xml.h>

using namespace dem;
using json = nlohmann::json;

void
BoundaryReader::read(const std::string& inputFileName, Box& container,
                     Box& grid, BoundaryPArray& boundaries) const
{
  std::ifstream ifs(inputFileName);
  if (!ifs) {
    std::cerr << "**ERROR**: Could not read boundary information from "
              << inputFileName << " in " << __FILE__ << ":" << __LINE__
              << std::endl;
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

  /*
  //std::cout << "container = " << container << "\n";
  //std::cout << "grid = " << grid << "\n";
  //std::cout << " Boundaries = \n";
  for (auto boundary : boundaries) {
    boundary->print(//std::cout);
    boundary->printContactInfo(//std::cout);
  }
  */
}

bool
BoundaryReader::readXML(const std::string& inputFileName, Box& container,
                        Box& grid, BoundaryPArray& boundaries) const
{
  // Read the input file
  zen::XmlDoc doc;
  try {
    //std::cout << "Input file name= " << inputFileName << "\n";
    doc = zen::load(inputFileName);
  } catch (const zen::XmlFileError& err) {
    std::cerr << "*ERROR** Could not read input file " << inputFileName << "\n";
    std::cerr << "    Error # = " << err.lastError << "\n";
    return false;
  } catch (const zen::XmlParsingError& err) {
    std::cerr << "*ERROR** Could not read input file " << inputFileName << "\n";
    std::cerr << "    Parse Error in line: " << err.row + 1
              << " col: " << err.col << "\n";
    return false;
  }

  // Check whether this is the right type of input file
  if (doc.root().getNameAs<std::string>() != "Ellip3D_input") {
    std::cerr << "*ERROR** Could not find tag <Ellip3D_input> in input file "
              << inputFileName << "\n";
    return false;
  }

  // Load the document into input proxy for easier element access
  zen::XmlIn ps(doc);

  // Read the title
  std::string title;
  if (!ps["Meta"]["title"](title)) {
    std::cerr << "*ERROR** Could not find boundary title in input file "
              << inputFileName << "\n";
    std::cerr << "  Add the <title> tag inside a <Meta> tag\n";
    return false;
  }
  //std::cout << "title = " << title << "\n";

  // Read the boundary information
  auto boundary_ps = ps["Boundary"];
  if (!boundary_ps) {
    std::cerr << "**ERROR** <Boundary> tag not found. \n";
    return false;
  }

  // Read the container dimensions
  std::string vecStr;
  if (!boundary_ps["containerMin"](vecStr)) {
    std::cerr
      << "**ERROR** Container min. position not found in boundary geometry\n";
    std::cerr << "  Add the <containerMin> [x, y, z] </containerMin> tag.";
    return false;
  }
  Vec boxMin = Vec::fromString(vecStr);

  if (!boundary_ps["containerMax"](vecStr)) {
    std::cerr
      << "**ERROR** Container max. position not found in boundary geometry\n";
    std::cerr << "  Add the <containerMax> [x, y, z] </containerMax> tag.";
    return false;
  }
  Vec boxMax = Vec::fromString(vecStr);

  container.set(boxMin.x(), boxMin.y(), boxMin.z(), boxMax.x(), boxMax.y(),
                boxMax.z());

  // compute grid assumed to be the same as container, change in
  // scatterParticle() if necessary.
  grid.set(boxMin.x(), boxMin.y(), boxMin.z(), boxMax.x(), boxMax.y(),
           boxMax.z());

  BoundaryType type;
  boundaries.clear();
  for (auto bound_ps = boundary_ps["boundary"]; bound_ps; bound_ps.next()) {
    std::string boundaryType;
    bound_ps.attribute("type", boundaryType);
    BoundaryId id;
    bound_ps.attribute("id", id);
    switch (getEnum(boundaryType)) {
      case PLANE:
        type = 1;
        boundaries.push_back(
          std::make_shared<PlaneBoundary>(id, type, bound_ps));
        break;
      case CYLINDER:
        type = 2;
        boundaries.push_back(
          std::make_shared<CylinderBoundary>(id, type, bound_ps));
        break;
      case NONE:
        break;
    }
  }

  /*
  //std::cout << "container = " << container << "\n";
  //std::cout << "grid = " << grid << "\n";
  //std::cout << " Boundaries = \n";
  for (auto boundary : boundaries) {
    boundary->print(//std::cout);
    boundary->printContactInfo(//std::cout);
  }
  */

  return true;
}

bool
BoundaryReader::readJSON(const std::string& inputFileName, Box& container,
                         Box& grid, BoundaryPArray& boundaries) const
{
  // Create an input ifstream
  std::ifstream ifs(inputFileName);
  if (!ifs) {
    std::cerr << "*ERROR** Could not read input file " << inputFileName << "\n";
    return false;
  }

  // Read the file stream into a string stream
  std::stringstream iss;
  iss << ifs.rdbuf();
  ifs.close();

  // Parse the input stream
  json doc;
  try {
    doc << iss;
    // //std::cout << std::setw(2) << doc << "\n";
  } catch (std::invalid_argument e) {
    std::cerr << "*ERROR** Could not parse input file " << inputFileName
              << "\n";
    std::cerr << "Please check for correctness using a linter.\n";
    return false;
  }

  // Check whether this is the right type of input file
  json ps;
  try {
    ps = doc["Ellip3D_input"];
  } catch (std::out_of_range e) {
    std::cerr << "*ERROR** Could not find key Ellip3D_input in input file "
              << inputFileName << "\n";
    return false;
  }

  // Read the title
  std::string title;
  try {
    auto meta_ps = ps["Meta"];
    auto title_ps = meta_ps.at("title");
    title = title_ps.get<std::string>();
  } catch (std::out_of_range e) {
    std::cerr << "*ERROR** Could not find boundary title in input file "
              << inputFileName << "\n";
    std::cerr << "  Add the \"title\" key inside a \"Meta\" object\n";
    return false;
  }
  //std::cout << "title = " << title << "\n";

  // Read the boundary information
  json boundary_ps;
  try {
    boundary_ps = ps["Boundary"];
  } catch (std::exception e) {
    std::cerr << "**ERROR** \"Boundary\" key not found. \n";
    return false;
  }

  // Read the container dimensions
  std::string vecStr;
  try {
    vecStr = boundary_ps["containerMin"].get<std::string>();
  } catch (std::exception e) {
    std::cerr
      << "**ERROR** Container min. position not found in boundary geometry\n";
    std::cerr << "  Add the containerMin: [x, y, z]  tag.";
    return false;
  }
  Vec boxMin = Vec::fromString(vecStr);

  try {
    vecStr = boundary_ps["containerMax"].get<std::string>();
  } catch (std::exception e) {
    std::cerr
      << "**ERROR** Container max. position not found in boundary geometry\n";
    std::cerr << "  Add the <containerMax> [x, y, z] </containerMax> tag.";
    return false;
  }
  Vec boxMax = Vec::fromString(vecStr);

  container.set(boxMin.x(), boxMin.y(), boxMin.z(), boxMax.x(), boxMax.y(),
                boxMax.z());

  // compute grid assumed to be the same as container, change in
  // scatterParticle() if necessary.
  grid.set(boxMin.x(), boxMin.y(), boxMin.z(), boxMax.x(), boxMax.y(),
           boxMax.z());

  BoundaryType type;
  boundaries.clear();
  try {
    auto bound_ps = boundary_ps["boundary"];
    // //std::cout << std::boolalpha << bound_ps.is_array() << "\n";
    // //std::cout << std::setw(2) << bound_ps << "\n";
    for (auto object : bound_ps) {
      std::string boundaryType = object["type"].get<std::string>();
      BoundaryId id = object["id"].get<BoundaryId>();
      switch (getEnum(boundaryType)) {
        case PLANE:
          type = 1;
          boundaries.push_back(
            std::make_shared<PlaneBoundary>(type, id, object));
          break;
        case CYLINDER:
          type = 2;
          boundaries.push_back(
            std::make_shared<CylinderBoundary>(type, id, object));
          break;
        case NONE:
          break;
      }
    }
  } catch (std::exception e) {
    std::cerr << "**ERROR** Boundaries not found in boundary geometry\n";
    std::cerr << "  Add the boundary key: value array tags.";
    return false;
  }

  return true;
}
