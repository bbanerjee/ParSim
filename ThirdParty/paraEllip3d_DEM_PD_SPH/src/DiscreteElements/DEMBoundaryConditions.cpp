#include <DiscreteElements/DEMBoundaryConditions.h>

using namespace dem;

bool 
DEMBoundaryConditions::read(const std::string& inputFileName)
{
  // Read the input file
  zen::XmlDoc doc;
  try {
    //std::cout << "Input file name= " << inputFileName << "\n";
    doc = zen::load(inputFileName);
  } catch (const zen::XmlFileError& err) {
    std::cerr << "*ERROR** Could not read boundary condition input file " 
              << inputFileName << "\n";
    std::cerr << "    Error # = " << err.lastError << "\n";
    return false;
  } catch (const zen::XmlParsingError& err) {
    std::cerr << "*ERROR** Could not read boundary condition input file " 
              << inputFileName << "\n";
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

  // Read the boundary condition information
  // **TODO** Add other types of BC
  auto periodic_bc_ps = ps["PeriodicParticleBC"];
  if (periodic_bc_ps) {
    d_domainBCType = BCUtils::DEM_DomainBCType::PERIODIC_PARTICLE;
    std::string bc_type;
    periodic_bc_ps.attribute("type", bc_type);
    if (bc_type == "deformation_gradient") {
      d_particleBCType = BCUtils::DEM_ParticleBCType::DEFORMATION_GRADIENT;
      d_particleDefGradBC.read(periodic_bc_ps);
    } else if (bc_type == "displacement") {
      d_particleBCType = BCUtils::DEM_ParticleBCType::DISPLACEMENT;
      d_particleDispBC.read(periodic_bc_ps);
    } else if (bc_type == "axisymmetric_strain") {
      d_particleBCType = BCUtils::DEM_ParticleBCType::AXISYMMETRIC_STRAIN;
      d_particleAxiStrainBC.read(periodic_bc_ps);
    }
  } else {
    d_domainBCType = BCUtils::DEM_DomainBCType::FIXED;
    d_particleBCType = BCUtils::DEM_ParticleBCType::NONE;
  }

  return true;
}
