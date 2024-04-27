#include <Core/Math/IntVec.h>
#include <Core/Math/Vec.h>
#include <Peridynamics/PeriParticle.h>
#include <Peridynamics/PeriElement.h>
#include <InputOutput/PeriParticleFileReader.h>
#include <fstream>
#include <sstream>
#include <type_traits>
#include <regex>

using namespace pd;
using dem::IntVec;
using dem::Vec;

void
PeriParticleFileReader::read(const std::string& fileName, 
                             PeriParticlePArray& particles,
                             PeriElementArray& connectivity) const
{
  // Try to open file
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::ostringstream out;
    out << "Could not open node input peridynamics volume mesh file " 
        << fileName << " for reading \n";
    std::cerr << out.str() << " in " << __FILE__ << " : " << __LINE__ << "\n";
    exit(-1);
  }
  file.close();

  // Check whether file is in Abaqus format
  bool isAbaqusFormat = checkAbaqusFileFormat(fileName);

  if (isAbaqusFormat) { 
    readPeriParticlesAbaqus(fileName, particles, connectivity);
  } else {
    readPeriParticlesText(fileName, particles, connectivity);
  }
}

bool 
PeriParticleFileReader::checkAbaqusFileFormat(const std::string& fileName) const
{
  // Read file
  std::ifstream file(fileName);
  std::string line;
  while (std::getline(file, line)) {

    // Ignore empty lines
    if (line.empty()) continue;

    // Erase white spaces from the line
    line.erase(remove(line.begin(),line.end(),' '),line.end());
    
    // Check whether the file contains a *Node directive
    std::transform (line.begin(), line.end(), line.begin(), ::tolower);
    if (line.compare("*node") == 0 ) {
      return true;
    }
  }
  return false;
}

/**
 *  Read particles directly from an input stream
 */
void
PeriParticleFileReader::readPeriParticlesText(const std::string& fileName,
                                              PeriParticlePArray& particles,
                                              PeriElementArray& connectivity) const
{
  std::ifstream ifs(fileName);

  int ndim = 0;
  int nPeriParticle = 0;
  int nele = 0;
  ifs >> ndim >> nPeriParticle >> nele;
  // read particle information, create and store PeriParticle objects into
  // periParticleVec
  for (int ip = 0; ip < nPeriParticle; ip++) {
    REAL tmp_x, tmp_y, tmp_z;
    ParticleID tmp_int;
    ifs >> tmp_int >> tmp_x >> tmp_y >> tmp_z;
    particles.push_back(
      std::make_shared<pd::PeriParticle>(tmp_int, tmp_x, tmp_y, tmp_z));
  }

  // read the connectivity information
  connectivity.resize(nele);
  for (int iel = 0; iel < nele; iel++) {
    ParticleID tmp_int;
    ifs >> tmp_int;
    for (ParticleID node = 0; node < 8; node++) {
      ifs >> connectivity[iel][node];
    }
  }

  ifs.close();
}

//--------------------------------------------------------------------------------
// Read input tetrahedral volume mesh file in Abaqus format
//--------------------------------------------------------------------------------
bool
PeriParticleFileReader::readPeriParticlesAbaqus(const std::string& fileName,
                                                PeriParticlePArray& particles,
                                                PeriElementArray& connectivity) const
{
  // Local Arrays
  std::vector<MeshNode> nodeArray;
  std::vector<VolumeElement> volElemArray;

  // Read file
  std::ifstream file(fileName);
  std::string line;
  bool node_flag = false;
  [[maybe_unused]] bool surf_elem_flag = false;
  [[maybe_unused]] bool line_elem_flag = false;
  bool vol_elem_flag = false;
  while (std::getline(file, line)) {

    // std::cout << "line = " << line << std::endl;

    // Ignore empty lines
    if (line.empty()) continue;

    // Erase white spaces from the line
    line.erase(remove(line.begin(),line.end(),' '),line.end());
    
    // Convert to lower case
    std::transform (line.begin(), line.end(), line.begin(), ::tolower);

    // Skip comment lines except *Node and *Element
    if (line[0] == '*') {
      node_flag = false;
      surf_elem_flag = false;
      vol_elem_flag = false;
      if (line.compare("*node") == 0) {
        node_flag = true;
      }
      std::string str("*element");
      if (line.compare(0, str.length(), str) == 0 ) {
        std::regex line_token(".*(line).*");
        std::regex surface_token(".*(surface).*");
        //std::cout << "line = " << line ;
        if (std::regex_search(line.begin(), line.end(), line_token)) {
          line_elem_flag = true;
        } else if (std::regex_search(line.begin(), line.end(), surface_token)) {
          //std::cout << " contains the string SURFACE";
          surf_elem_flag = true;
        } else {
          //std::cout << " does not contain the string SURFACE: Volume element";
          vol_elem_flag = true;
        }
        //std::cout << std::endl;
      }
      std::regex end_token(".*(end).*");
      if (std::regex_match(line.begin(), line.end(), end_token)) {
        //std::cout << "Line contains *End" << std::endl;
        node_flag = false;
        surf_elem_flag = false;
        vol_elem_flag = false;
        line_elem_flag = false;
      }
      continue;
    }

    // Read the nodal coordinates
    if (node_flag) {
      readAbaqusMeshNode(line, nodeArray);
    }

    // Read the volume element connectivity
    if (vol_elem_flag) {
      readAbaqusMeshVolumeElement(line, volElemArray);
    }

  }
  file.close();

  //std::cout << "Done reading " << std::endl;

  // Sort volume elements in ascending order by ID
  std::sort(volElemArray.begin(), volElemArray.end(), 
            [](const VolumeElement& a, const VolumeElement& b) {
              return b.id_ > a.id_;
            });

  // Print volume elements
  /*
  std::cout << "Volume elements : " << std::endl;
  for (auto iter = volElemArray.begin(); iter != volElemArray.end(); iter++) {
    VolumeElement volElem = *iter;
    std::cout << "\t ID = " << volElem.id_
              << ", nodes = (" << volElem.node1_ << ", " << volElem.node2_ 
              << ", " << volElem.node3_ << ", " << volElem.node4_ << ")" 
              << std::endl;
  }

  std::cout << "Nodes : " << std::endl;
  for (auto iter = nodeArray.begin(); iter != nodeArray.end(); iter++) {
    MeshNode node = *iter;
    std::cout << "\t ID = " << node.id_
              << ", pos = (" << node.x_ << ", " << node.y_ << ", " << node.z_
              << ")" << std::endl;
  }
  */

  // Now push the coordinates 
  for (auto& node : nodeArray) {
    particles.push_back(
      std::make_shared<pd::PeriParticle>(node.id_, node.x_, node.y_, node.z_));
  }
  // read the connectivity information
  connectivity.resize(volElemArray.size());
  int elemID = 0;
  for (auto& element : volElemArray) {
    connectivity[elemID][0] = element.node1_;
    connectivity[elemID][1] = element.node2_;
    connectivity[elemID][2] = element.node3_;
    connectivity[elemID][3] = element.node4_;
    elemID++;
  }

  file.close();

  return true;
}

void
PeriParticleFileReader::readAbaqusMeshNode(const std::string& inputLine,
                                           std::vector<MeshNode>& nodes) const
{
  // Tokenize the string
  std::string data_piece;
  std::istringstream data_stream(inputLine);
  std::vector<std::string> data;
  while (std::getline(data_stream, data_piece, ',')) {
    data.push_back(data_piece);
  }

  if (data.size() < 4) {
    std::ostringstream out;
    out << "Could not read nodal coordinates from " 
        << inputLine << std::endl;
    std::cerr << out.str() << " in " << __FILE__ << " : " << __LINE__ << "\n";
    return;
  }

  auto iter = data.begin(); 
  ParticleID node_id = std::stoi(*iter); ++iter;
  double xcoord = std::stod(*iter); ++iter;
  double ycoord = std::stod(*iter); ++iter;
  double zcoord = std::stod(*iter);

  // Save the nodal coordinates
  nodes.emplace_back(MeshNode(node_id, xcoord, ycoord, zcoord));
}

void
PeriParticleFileReader::readAbaqusMeshVolumeElement(const std::string& inputLine,
                                                    std::vector<VolumeElement>& elements) const
{
  // Tokenize the string
  std::string data_piece;
  std::istringstream data_stream(inputLine);
  std::vector<std::string> data;
  while (std::getline(data_stream, data_piece, ',')) {
    data.push_back(data_piece);
  }

  if (data.size() != 5) {
    std::ostringstream out;
    out << "Could not read volume element connectivity from input line: " 
        << inputLine << std::endl;
    std::cerr << out.str() << " in " << __FILE__ << " : " << __LINE__ << "\n";
    return;
  }

  // Read the element id
  auto iter = data.begin(); 
  ElementID element_id = std::stoi(*iter); ++iter;

  // Read the element node ids
  std::vector<ParticleID> node_list;
  for (; iter != data.end(); ++iter) {
    ParticleID node_id = std::stoi(*iter);
    node_list.emplace_back(node_id);
  }
  if (node_list.empty()) {
    std::ostringstream out;
    out << "Could not find nodes in volume element input data stream.";
    std::cerr << out.str() << " in " << __FILE__ << " : " << __LINE__ << "\n";
    return;
  }

  // Save the data
  elements.emplace_back(VolumeElement(element_id, node_list));
}
