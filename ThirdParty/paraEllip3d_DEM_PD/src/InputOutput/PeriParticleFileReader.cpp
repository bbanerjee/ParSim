#include "zlib.h"
#include <Core/Math/IntVec.h>
#include <Core/Math/Vec.h>
#include <Peridynamics/PeriParticle.h>
#include <Peridynamics/PeriElement.h>
#include <InputOutput/PeriParticleFileReader.h>
#include <InputOutput/cppcodec/cppcodec/base64_default_rfc4648.hpp>
#include <fstream>
#include <sstream>
#include <type_traits>

using namespace pd;
using dem::IntVec;
using dem::Vec;

void
PeriParticleFileReader::read(const std::string& fileName, 
                             PeriParticlePArray& particles,
                             PeriElements& connectivity) const
{
  // Check whether file is XML
  std::ifstream ifs(fileName);
  char firstChar;
  ifs >> firstChar;
  ifs.close();
  if (firstChar == '<') { // XML
    readPeriParticlesXML(fileName, particles, connectivity);
  } else {
    readPeriParticlesText(fileName, particles, connectivity);
  }
}

/**
 *  Read particles directly from an input stream
 */
void
PeriParticleFileReader::readPeriParticlesText(const std::string& inputFileName,
                                              PeriParticlePArray& particles,
                                              PeriElements& connectivity) const
{
  std::ifstream ifs(inputFileName);
  if (!ifs) {
    std::cerr << "**ERROR**: Could not read input particle file  "
              << inputFileName << " in " << __FILE__ << ":" << __LINE__
              << std::endl;
    exit(-1);
  }

  int ndim = 0;
  int nPeriParticle = 0;
  int nele = 0;
  ifs >> ndim >> nPeriParticle >> nele;
  // read particle information, create and store PeriParticle objects into
  // periParticleVec
  for (int ip = 0; ip < nPeriParticle; ip++) {
    REAL tmp_x, tmp_y, tmp_z;
    int tmp_int;
    ifs >> tmp_int >> tmp_x >> tmp_y >> tmp_z;
    particles.push_back(
      std::make_shared<pd::PeriParticle>(tmp_x, tmp_y, tmp_z));
  }

  // read the connectivity information
  for (int iel = 0; iel < nele; iel++) {
    int tmp_int;
    ifs >> tmp_int;
    for (int node = 0; node < 8; node++) {
      ifs >> connectivity[iel][node];
    }
  }

  ifs.close();
}

bool
PeriParticleFileReader::readPeriParticlesXML(const std::string& inputFileName,
                                             PeriParticlePArray& particles,
                                             PeriElements& connectivity) const
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
      bool success = readPeriParticleValues<size_t>(particle_ps, "id", particleType,
                                                particleIDs);

      // Get the particle radii
      std::vector<Vec> particleRadii;
      success = readPeriParticleValues<Vec>(particle_ps, "radii", particleType,
                                        particleRadii);

      // Get the particle axle_a, axle_b, axle_c
      std::vector<Vec> particleAxleA, particleAxleB, particleAxleC;
      success = readPeriParticleValues<Vec>(particle_ps, "axle_a", particleType,
                                        particleAxleA);
      success = readPeriParticleValues<Vec>(particle_ps, "axle_b", particleType,
                                        particleAxleB);
      success = readPeriParticleValues<Vec>(particle_ps, "axle_c", particleType,
                                        particleAxleC);

      // Get the particle position, velocity, omega
      std::vector<Vec> particlePos, particleVel, particleRot;
      success = readPeriParticleValues<Vec>(particle_ps, "position", particleType,
                                        particlePos);
      success = readPeriParticleValues<Vec>(particle_ps, "velocity", particleType,
                                        particleVel);
      success = readPeriParticleValues<Vec>(particle_ps, "omega", particleType,
                                        particleRot);

      // Get the particle force and moment
      std::vector<Vec> particleForce, particleMoment;
      success = readPeriParticleValues<Vec>(particle_ps, "force", particleType,
                                        particleForce);
      success = readPeriParticleValues<Vec>(particle_ps, "moment", particleType,
                                        particleMoment);

      if (!success) {
        std::cerr << "Read failed\n";
        return false;
      }

      // **TODO** Assign type num using particleType
      std::size_t particleTypeNum = 0;
      for (std::size_t ii = 0; ii < numParticles; ++ii) {
        PeriParticleP pt = std::make_shared<PeriParticle>(
          particlePos[ii].x(), particlePos[ii].y(), particlePos[ii].z());

        particles.push_back(pt);
      }
    } else {
      std::cerr << "Logic for ASCII text files not implemented yet\n";
      return false;
    }

  } // end of loop over particles

  return true;
}

template <typename T>
bool
PeriParticleFileReader::readPeriParticleValues(zen::XmlIn& ps, const std::string& name,
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
PeriParticleFileReader::decodeAndUncompress(const std::string& inputStr,
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
PeriParticleFileReader::convert<int>(const std::string& str) const
{
  return std::stoi(str);
}

template <>
size_t
PeriParticleFileReader::convert<size_t>(const std::string& str) const
{
  return std::stoul(str);
}

template <>
double
PeriParticleFileReader::convert<double>(const std::string& str) const
{
  return std::stod(str);
}

template <>
Vec
PeriParticleFileReader::convert<Vec>(const std::string& str) const
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  return Vec(std::stod(split[0]), std::stod(split[1]), std::stod(split[2]));
}

template <>
IntVec
PeriParticleFileReader::convert<IntVec>(const std::string& str) const
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = { std::istream_iterator<std::string>{ iss },
                                     std::istream_iterator<std::string>{} };
  return IntVec(std::stoi(split[0]), std::stoi(split[1]), std::stoi(split[2]));
}

template <>
std::vector<double>
PeriParticleFileReader::convert<std::vector<double>>(const std::string& str) const
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
PeriParticleFileReader::convertStrArray(const std::string& str) const
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
