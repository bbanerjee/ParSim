#include "zlib.h"
#include <Core/Math/IntVec.h>
#include <Core/Math/Vec.h>
#include <DiscreteElements/Particle.h>
#include <InputOutput/ParticleFileReader.h>
#include <InputOutput/cppcodec/cppcodec/base64_default_rfc4648.hpp>
#include <fstream>
#include <sstream>
#include <type_traits>

using namespace dem;

void
ParticleFileReader::read(const std::string& fileName, const REAL& youngModulus,
                         const REAL& poissonRatio, bool doInitialize,
                         ParticlePArray& particles, Gradation& gradation)
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
ParticleFileReader::readParticlesText(const std::string& inputParticle,
                                      ParticlePArray& allParticleVec,
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

  allParticleVec.clear();

  std::size_t id, type;
  REAL a, b, c, px, py, pz, dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz;
  REAL vx, vy, vz, omx, omy, omz, fx, fy, fz, mx, my, mz;
  // REAL young = util::getParam<REAL>("young");
  // REAL poisson = util::getParam<REAL>("poisson");
  for (std::size_t i = 0; i < particleNum; ++i) {
    ifs >> id >> type >> a >> b >> c >> px >> py >> pz >> dax >> day >> daz >>
      dbx >> dby >> dbz >> dcx >> dcy >> dcz >> vx >> vy >> vz >> omx >> omy >>
      omz >> fx >> fy >> fz >> mx >> my >> mz;

    ParticleP pt = std::make_shared<Particle>(
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
      pt->setPrevOmga(Vec(omx, omy, omz));
      pt->setCurrOmga(Vec(omx, omy, omz));
      pt->setForce(Vec(fx, fy, fz));  // initial force
      pt->setMoment(Vec(mx, my, mz)); // initial moment
    }

    allParticleVec.push_back(pt);
  }

  std::size_t sieveNum;
  ifs >> sieveNum;
  std::vector<REAL> percent(sieveNum), size(sieveNum);
  for (auto i = 0u; i < sieveNum; ++i)
    ifs >> percent[i] >> size[i];
  REAL ratio_ba, ratio_ca;
  ifs >> ratio_ba >> ratio_ca;
  gradation.set(sieveNum, percent, size, ratio_ba, ratio_ca);

  ifs.close();
}

bool
ParticleFileReader::readParticlesXML(const std::string& inputFileName,
                                     ParticlePArray& allParticleVec,
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
        ParticleP pt = std::make_shared<Particle>(
          particleIDs[ii], particleTypeNum, particleRadii[ii], particlePos[ii],
          particleAxleA[ii], particleAxleB[ii], particleAxleC[ii],
          d_youngModulus, d_poissonRatio);

        // optional settings for a particle's initial status
        // std::cout << "doInitialize = " << std::boolalpha << d_doInitialize <<
        // std::endl;
        if (d_doInitialize) {
          pt->setPrevVeloc(particleVel[ii]);
          pt->setCurrVeloc(particleVel[ii]);
          pt->setPrevOmga(particleRot[ii]);
          pt->setCurrOmga(particleRot[ii]);
          pt->setForce(particleForce[ii]);   // initial force
          pt->setMoment(particleMoment[ii]); // initial moment
        }

        allParticleVec.push_back(pt);
      }
    } else {
      std::cerr << "Logic for ASCII text files not implemented yet\n";
      return false;
    }

  } // end of loop over particles

  // Sieve data (assume ASCII)
  // **TODO** Add validity checks
  auto sieve_ps = ps["Sieves"];
  if (sieve_ps) {
    std::size_t numSieves;
    sieve_ps.attribute("number", numSieves);

    std::string percentPassingStr;
    sieve_ps["percent_passing"](percentPassingStr);
    std::vector<REAL> percentPassing = convertStrArray<REAL>(percentPassingStr);

    std::string sizeStr;
    sieve_ps["size"](sizeStr);
    std::vector<REAL> size = convertStrArray<REAL>(sizeStr);

    REAL ratio_ba, ratio_ca;
    sieve_ps["sieve_ratio"]["ratio_ba"](ratio_ba);
    sieve_ps["sieve_ratio"]["ratio_ca"](ratio_ca);

    gradation.set(numSieves, percentPassing, size, ratio_ba, ratio_ca);
  }

  return true;
}

namespace dem {

template <typename T>
bool
ParticleFileReader::readParticleValues(zen::XmlIn& ps, const std::string& name,
                                       const std::string& particleType,
                                       std::vector<T>& output) const
{
  // Get the particle values
  auto prop_ps = ps[name];
  if (!prop_ps) {
    std::cerr << "**ERROR** Particle " << name << " not found "
              << " particle type: " << particleType << "\n";
    return false;
  }
  std::string particleDataStr;
  prop_ps(particleDataStr);
  int numComp = 0;
  prop_ps.attribute("numComponents", numComp);
  bool success = decodeAndUncompress<T>(particleDataStr, numComp, output);
  if (!success) {
    std::cerr << "**ERROR** Could not decode and uncompress particle " << name
              << "\n";
    return false;
  }
  return true;
}

template <typename T>
bool
ParticleFileReader::decodeAndUncompress(const std::string& inputStr,
                                        const int& numComponents,
                                        std::vector<T>& output) const
{
  // Decode from base64
  std::vector<std::uint8_t> decoded = base64::decode(inputStr);

  // Uncompress from gzip
  std::vector<std::uint8_t> uncompressed;
  z_stream stream;

  // Allocate inflate state
  stream.zalloc = Z_NULL;
  stream.zfree = Z_NULL;
  stream.opaque = Z_NULL;
  stream.avail_in = 0;
  stream.next_in = Z_NULL;
  int err = inflateInit(&stream);
  if (err != Z_OK) {
    std::cerr << "inflateInit"
              << " error: " << err << std::endl;
    return false;
  }

  // Uncompress until stream ends
  stream.avail_in = decoded.size();
  stream.next_in = &decoded[0];
  do {
    do {
      std::vector<std::uint8_t> out(decoded.size());
      stream.avail_out = out.size();
      stream.next_out = &out[0];

      err = inflate(&stream, Z_SYNC_FLUSH);
      // auto have = decoded.size() - stream.avail_out;
      uncompressed.insert(std::end(uncompressed), std::begin(out),
                          std::end(out));

    } while (stream.avail_out == 0);
  } while (err != Z_STREAM_END);

  if (inflateEnd(&stream) != Z_OK) {
    std::cerr << "inflateEnd"
              << " error: " << err << std::endl;
    return false;
  }

  // Split the uncompressed string into a vector of tokens
  // (Assume that data are space separated)
  // (See: https://stackoverflow.com/questions/236129/split-a-string-in-c)
  std::istringstream iss(std::string(uncompressed.begin(), uncompressed.end()));
  std::vector<std::string> outputStr = { std::istream_iterator<std::string>{
                                           iss },
                                         std::istream_iterator<std::string>{} };
  for (auto iter = outputStr.begin(); iter != outputStr.end();
       iter += numComponents) {
    std::string str = *iter;
    for (int ii = 1; ii < numComponents; ii++) {
      str += " ";
      str += *(iter + ii);
    }
    output.push_back(convert<T>(str));
  }

  return true;
}

template <>
int
ParticleFileReader::convert<int>(const std::string& str) const
{
  return std::stoi(str);
}

template <>
size_t
ParticleFileReader::convert<size_t>(const std::string& str) const
{
  return std::stoul(str);
}

template <>
double
ParticleFileReader::convert<double>(const std::string& str) const
{
  return std::stod(str);
}

template <>
Vec
ParticleFileReader::convert<Vec>(const std::string& str) const
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  return Vec(std::stod(split[0]), std::stod(split[1]), std::stod(split[2]));
}

template <>
IntVec
ParticleFileReader::convert<IntVec>(const std::string& str) const
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  return IntVec(std::stoi(split[0]), std::stoi(split[1]), std::stoi(split[2]));
}

template <>
std::vector<double>
ParticleFileReader::convert<std::vector<double>>(const std::string& str) const
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  std::vector<double> vec;
  for (auto str : split) {
    vec.push_back(std::stod(str));
  }
  return vec;
}

template <typename T>
std::vector<T>
ParticleFileReader::convertStrArray(const std::string& str) const
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  std::vector<T> vec;

  if (std::is_same<T, REAL>::value) {
    for (auto str : split) {
      vec.push_back(std::stod(str));
    }
  } else if (std::is_same<T, size_t>::value) {
    for (auto str : split) {
      vec.push_back(std::stoul(str));
    }
  } else {
    std::cerr << "**ERROR** Conversion of string array is allowed only for"
              << " numeric types\n";
  }
  return vec;
}
} // end namespace dem