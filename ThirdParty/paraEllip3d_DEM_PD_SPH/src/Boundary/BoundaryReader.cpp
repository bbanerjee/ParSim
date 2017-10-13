#include <Boundary/BoundaryReader.h>
#include <Boundary/BoundaryContainers.h>
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
                     Box& patchBox, BoundaryPArray& boundaries) const
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
  // compute patchGrid assumed to be the same as container, change in
  // scatterParticles() if necessary.
  patchBox.set(x1, y1, z1, x2, y2, z2);

  boundaries.clear();
  std::size_t boundaryNum;
  ifs >> boundaryNum;
  for (std::size_t i = 0; i < boundaryNum; ++i) {
    std::size_t type;
    ifs >> type;
    Boundary::BoundaryType boundaryType = Boundary::getBoundaryType(type);
    if (type == 1) // plane boundary
      boundaries.push_back(std::make_shared<PlaneBoundary>(boundaryType, ifs));
    else if (type == 2) // cylindrical boundary
      boundaries.push_back(std::make_shared<CylinderBoundary>(boundaryType, ifs));
  }

  ifs.close();

  /*
  std::cout << "container = " << container << "\n";
  std::cout << "patchBox = " << patchBox << "\n";
  //std::cout << " Boundaries = \n";
  for (auto boundary : boundaries) {
    boundary->print(//std::cout);
    boundary->printContactInfo(//std::cout);
  }
  */
}

bool
BoundaryReader::readXML(const std::string& inputFileName, Box& container,
                        Box& patchBox, BoundaryPArray& boundaries) const
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

  // compute patchGrid assumed to be the same as container, change in
  // scatterParticles() if necessary.
  patchBox.set(boxMin.x(), boxMin.y(), boxMin.z(), boxMax.x(), boxMax.y(),
           boxMax.z());

  boundaries.clear();
  for (auto bound_ps = boundary_ps["boundary"]; bound_ps; bound_ps.next()) {

    std::string boundaryType;
    bound_ps.attribute("type", boundaryType);
    Boundary::BoundaryType type = Boundary::getBoundaryType(boundaryType);

    std::string boundaryID;
    bound_ps.attribute("id", boundaryID);
    Boundary::BoundaryID id = Boundary::getBoundaryID(boundaryID);

    switch (type) {
      case Boundary::BoundaryType::PLANE:
        boundaries.push_back(
          std::make_shared<PlaneBoundary>(type, id, bound_ps));
        break;
      case Boundary::BoundaryType::CYLINDER:
        boundaries.push_back(
          std::make_shared<CylinderBoundary>(type, id, bound_ps));
        break;
      case Boundary::BoundaryType::NONE:
        break;
    }
  }

  /*
  std::cout << "container = " << container << "\n";
  std::cout << "patchBox = " << patchBox << "\n";
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
                         Box& patchBox, BoundaryPArray& boundaries) const
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

  // compute patchGrid assumed to be the same as container, change in
  // scatterParticles() if necessary.
  patchBox.set(boxMin.x(), boxMin.y(), boxMin.z(), boxMax.x(), boxMax.y(),
           boxMax.z());

  boundaries.clear();
  try {
    auto bound_ps = boundary_ps["boundary"];
    // //std::cout << std::boolalpha << bound_ps.is_array() << "\n";
    // //std::cout << std::setw(2) << bound_ps << "\n";
    for (auto object : bound_ps) {
      std::string boundaryType = object["type"].get<std::string>();
      std::string boundaryID = object["id"].get<std::string>();
      Boundary::BoundaryType type = Boundary::getBoundaryType(boundaryType);
      Boundary::BoundaryID id = Boundary::getBoundaryID(boundaryID);
      switch (type) {
        case Boundary::BoundaryType::PLANE:
          boundaries.push_back(
            std::make_shared<PlaneBoundary>(type, id, object));
          break;
        case Boundary::BoundaryType::CYLINDER:
          boundaries.push_back(
            std::make_shared<CylinderBoundary>(type, id, object));
          break;
        case Boundary::BoundaryType::NONE:
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
