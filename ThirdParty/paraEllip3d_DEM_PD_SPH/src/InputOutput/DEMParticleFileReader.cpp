#include <InputOutput/DEMParticleFileReader.h>
#include <InputOutput/IOUtils.h>
#include <Core/Math/IntVec.h>
#include <Core/Math/Vec.h>
#include <DiscreteElements/DEMParticle.h>
#include <fstream>
#include <sstream>
#include <type_traits>

using namespace dem;

void
DEMParticleFileReader::read(const std::string& fileName, const REAL& youngModulus,
                         const REAL& poissonRatio, bool doInitialize,
                         DEMParticlePArray& particles, Gradation& gradation)
{
  d_youngModulus = youngModulus;
  d_poissonRatio = poissonRatio;
  d_doInitialize = doInitialize;

  // Check whether file is XML
  std::ifstream ifs(fileName);
  char firstChar;
  ifs >> firstChar;
  ifs.close();
  if (firstChar == '<') { // XML
    readParticlesXML(fileName, particles, gradation);
  } else {
    readParticlesText(fileName, particles, gradation);
  }
}

/**
 *  Read particles directly from an input stream
 */
void
DEMParticleFileReader::readParticlesText(const std::string& inputParticle,
                                      DEMParticlePArray& allDEMParticleVec,
                                      Gradation& gradation) const
{
  std::ifstream ifs(inputParticle);
  if (!ifs) {
    std::cerr << "**ERROR**: Could not read input particle file  "
              << inputParticle << " in " << __FILE__ << ":" << __LINE__
              << std::endl;
    exit(-1);
  }

  std::size_t particleNum;
  ifs >> particleNum;

  std::string str;
  ifs >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >>
    str >> str >> str >> str >> str >> str >> str >> str >> str >> str >> str >>
    str >> str >> str >> str >> str >> str >> str >> str;

  allDEMParticleVec.clear();

  std::size_t id, type;
  REAL a, b, c, px, py, pz, dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz;
  REAL vx, vy, vz, omx, omy, omz, fx, fy, fz, mx, my, mz;
  // REAL young = util::getParam<REAL>("young");
  // REAL poisson = util::getParam<REAL>("poisson");
  for (std::size_t i = 0; i < particleNum; ++i) {
    ifs >> id >> type >> a >> b >> c >> px >> py >> pz >> dax >> day >> daz >>
      dbx >> dby >> dbz >> dcx >> dcy >> dcz >> vx >> vy >> vz >> omx >> omy >>
      omz >> fx >> fy >> fz >> mx >> my >> mz;

    DEMParticleP pt = std::make_shared<DEMParticle>(
      id, type, Vec(a, b, c), Vec(px, py, pz), Vec(dax, day, daz),
      Vec(dbx, dby, dbz), Vec(dcx, dcy, dcz), d_youngModulus, d_poissonRatio);

    // optional settings for a particle's initial status
    // if ((static_cast<std::size_t>(
    //      "toInitParticle"])) == 1) {
    // //std::cout << "doInitialize = " << std::boolalpha << d_doInitialize <<
    // std::endl;
    if (d_doInitialize) {
      pt->setPrevVeloc(Vec(vx, vy, vz));
      pt->setCurrVeloc(Vec(vx, vy, vz));
      pt->setPrevOmega(Vec(omx, omy, omz));
      pt->setCurrOmega(Vec(omx, omy, omz));
      pt->setForce(Vec(fx, fy, fz));  // initial force
      pt->setMoment(Vec(mx, my, mz)); // initial moment
    }

    allDEMParticleVec.push_back(pt);
  }

  // Save the gradation data
  gradation.initializeFromCSVFile(ifs);

  ifs.close();
}

bool
DEMParticleFileReader::readParticlesXML(const std::string& inputFileName,
                                     DEMParticlePArray& allDEMParticleVec,
                                     Gradation& gradation) const
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

  // Loop through the particle types in the input file
  for (auto particle_ps = ps["Particles"]; particle_ps; particle_ps.next()) {

    // Get the attributes of the particles
    std::size_t numParticles = 0;
    particle_ps.attribute("number", numParticles);
    std::string particleType = "sphere";
    particle_ps.attribute("type", particleType);
    std::string compression = "none";
    particle_ps.attribute("compression", compression);
    std::string encoding = "none";
    particle_ps.attribute("encoding", encoding);

    if (encoding == "base64" && compression == "gzip") {

      // Get the particle ids
      std::vector<size_t> particleIDs;
      bool success = readParticleValues<size_t>(particle_ps, "id", particleType,
                                                particleIDs);

      // Get the particle radii
      std::vector<Vec> particleRadii;
      success = readParticleValues<Vec>(particle_ps, "radii", particleType,
                                        particleRadii);

      // Get the particle axle_a, axle_b, axle_c
      std::vector<Vec> particleAxleA, particleAxleB, particleAxleC;
      success = readParticleValues<Vec>(particle_ps, "axle_a", particleType,
                                        particleAxleA);
      success = readParticleValues<Vec>(particle_ps, "axle_b", particleType,
                                        particleAxleB);
      success = readParticleValues<Vec>(particle_ps, "axle_c", particleType,
                                        particleAxleC);

      // Get the particle position, velocity, omega
      std::vector<Vec> particlePos, particleVel, particleRot;
      success = readParticleValues<Vec>(particle_ps, "position", particleType,
                                        particlePos);
      success = readParticleValues<Vec>(particle_ps, "velocity", particleType,
                                        particleVel);
      success = readParticleValues<Vec>(particle_ps, "omega", particleType,
                                        particleRot);

      // Get the particle force and moment
      std::vector<Vec> particleForce, particleMoment;
      success = readParticleValues<Vec>(particle_ps, "force", particleType,
                                        particleForce);
      success = readParticleValues<Vec>(particle_ps, "moment", particleType,
                                        particleMoment);

      if (!success) {
        std::cerr << "Read failed\n";
        return false;
      }

      // **TODO** Assign type num using particleType
      std::size_t particleTypeNum = 0;
      for (std::size_t ii = 0; ii < numParticles; ++ii) {
        DEMParticleP pt = std::make_shared<DEMParticle>(
          particleIDs[ii], particleTypeNum, particleRadii[ii], particlePos[ii],
          particleAxleA[ii], particleAxleB[ii], particleAxleC[ii],
          d_youngModulus, d_poissonRatio);

        // optional settings for a particle's initial status
        // std::cout << "doInitialize = " << std::boolalpha << d_doInitialize <<
        // std::endl;
        if (d_doInitialize) {
          pt->setPrevVeloc(particleVel[ii]);
          pt->setCurrVeloc(particleVel[ii]);
          pt->setPrevOmega(particleRot[ii]);
          pt->setCurrOmega(particleRot[ii]);
          pt->setForce(particleForce[ii]);   // initial force
          pt->setMoment(particleMoment[ii]); // initial moment
        }

        allDEMParticleVec.push_back(pt);
      }
    } else {
      std::cerr << "Logic for ASCII text files not implemented yet\n";
      return false;
    }

  } // end of loop over particles

  // Save the gradation data
  gradation.initializeFromXMLFile(ps);

  return true;
}

namespace dem {

template <typename T>
bool
DEMParticleFileReader::readParticleValues(zen::XmlIn& ps, const std::string& name,
                                       const std::string& particleType,
                                       std::vector<T>& output) const
{
  // Get the particle values
  auto prop_ps = ps[name];
  if (!prop_ps) {
    std::cerr << "**ERROR** DEMParticle " << name << " not found "
              << " particle type: " << particleType << "\n";
    return false;
  }
  std::string particleDataStr;
  prop_ps(particleDataStr);
  int numComp = 0;
  prop_ps.attribute("numComponents", numComp);
  bool success = 
    Ellip3D::Util::decodeAndUncompress<T>(particleDataStr, numComp, output);
  if (!success) {
    std::cerr << "**ERROR** Could not decode and uncompress particle " << name
              << "\n";
    return false;
  }
  return true;
}

} // end namespace dem