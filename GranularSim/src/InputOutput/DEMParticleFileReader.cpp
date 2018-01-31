#include <InputOutput/DEMParticleFileReader.h>
#include <InputOutput/IOUtils.h>
#include <Core/Math/IntVec.h>
#include <Core/Math/Vec.h>
#include <DiscreteElements/DEMParticle.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
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

  // Compute total mass and save in each particle
  updateTotalMass(particles);
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
      id, 
      DEMParticle::DEMParticleShape::ELLIPSOID, 
      static_cast<DEMParticle::DEMParticleType>(type),
      Vec(a, b, c), Vec(px, py, pz), Vec(dax, day, daz),
      Vec(dbx, dby, dbz), Vec(dcx, dcy, dcz), d_youngModulus, d_poissonRatio);

    // optional settings for a particle's initial status
    // if ((static_cast<std::size_t>(
    //      "toInitParticle"])) == 1) {
    // //std::cout << "doInitialize = " << std::boolalpha << d_doInitialize <<
    // std::endl;
    if (d_doInitialize) {
      pt->setPreviousVelocity(Vec(vx, vy, vz));
      pt->setCurrentVelocity(Vec(vx, vy, vz));
      pt->setPreviousAngularVelocity(Vec(omx, omy, omz));
      pt->setCurrentAngularVelocity(Vec(omx, omy, omz));
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
    std::string particleShapeStr = "sphere";
    particle_ps.attribute("shape", particleShapeStr);
    std::string compression = "none";
    particle_ps.attribute("compression", compression);
    std::string encoding = "none";
    particle_ps.attribute("encoding", encoding);

    if (encoding == "base64" && compression == "gzip") {

      // Get the particle ids
      std::vector<size_t> particleIDs;
      bool success = readParticleValues<size_t>(particle_ps, "id", particleShapeStr,
                                                particleIDs);

      // Get the particle types
      std::vector<size_t> particleTypes;
      success = readParticleValues<size_t>(particle_ps, "type", particleShapeStr,
                                           particleTypes);

      // Get the particle radii
      std::vector<Vec> particleRadii;
      success = readParticleValues<Vec>(particle_ps, "radii", particleShapeStr,
                                        particleRadii);

      // Get the particle axle_a, axle_b, axle_c
      std::vector<Vec> particleAxleA, particleAxleB, particleAxleC;
      success = readParticleValues<Vec>(particle_ps, "axle_a", particleShapeStr,
                                        particleAxleA);
      success = readParticleValues<Vec>(particle_ps, "axle_b", particleShapeStr,
                                        particleAxleB);
      success = readParticleValues<Vec>(particle_ps, "axle_c", particleShapeStr,
                                        particleAxleC);

      // Get the particle position, velocity, omega
      std::vector<Vec> particlePos, particleVel, particleRot;
      success = readParticleValues<Vec>(particle_ps, "position", particleShapeStr,
                                        particlePos);
      success = readParticleValues<Vec>(particle_ps, "velocity", particleShapeStr,
                                        particleVel);
      success = readParticleValues<Vec>(particle_ps, "omega", particleShapeStr,
                                        particleRot);

      // Get the particle force and moment
      std::vector<Vec> particleForce, particleMoment;
      success = readParticleValues<Vec>(particle_ps, "force", particleShapeStr,
                                        particleForce);
      success = readParticleValues<Vec>(particle_ps, "moment", particleShapeStr,
                                        particleMoment);

      if (!success) {
        std::cerr << "Read failed\n";
        return false;
      }

      auto particleShape = DEMParticle::getDEMParticleShape(particleShapeStr);
      for (std::size_t ii = 0; ii < numParticles; ++ii) {
        DEMParticleP pt = std::make_shared<DEMParticle>(
          particleIDs[ii], 
          particleShape, 
          static_cast<DEMParticle::DEMParticleType>(particleTypes[ii]), 
          particleRadii[ii], particlePos[ii],
          particleAxleA[ii], particleAxleB[ii], particleAxleC[ii],
          d_youngModulus, d_poissonRatio);

        // optional settings for a particle's initial status
        // std::cout << "doInitialize = " << std::boolalpha << d_doInitialize <<
        // std::endl;
        if (d_doInitialize) {
          pt->setPreviousVelocity(particleVel[ii]);
          pt->setCurrentVelocity(particleVel[ii]);
          pt->setPreviousAngularVelocity(particleRot[ii]);
          pt->setCurrentAngularVelocity(particleRot[ii]);
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

void
DEMParticleFileReader::updateTotalMass(DEMParticlePArray& particles)
{
  double totalMass = 0.0;
  for (const auto& particle : particles) {
    totalMass += particle->mass();
  }

  for (auto& particle : particles) {
    particle->setTotalMass(totalMass);
  }
}

using vtkUnstructuredGridReaderP = vtkSmartPointer<vtkXMLUnstructuredGridReader>;
void 
DEMParticleFileReader::readVTK(const std::string& filename, 
                               REAL youngModulus, REAL poissonRatio,
                               DEMParticlePArray& particles)
{
  // Read all the data from the file
  vtkUnstructuredGridReaderP reader = vtkUnstructuredGridReaderP::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  // Get a pointer to the data
  vtkUnstructuredGrid* data = reader->GetOutput();

  // Get the point coordinates
  std::vector<Vec> positions;
  auto numPts = data->GetNumberOfPoints();
  for (int ii = 0; ii < numPts; ++ii) {
    double* pt = data->GetPoint(ii);
    positions.push_back(Vec(pt[0], pt[1], pt[2]));
    //std::cout << "(" << pt[0] << "," << pt[1] << "," << pt[2] << ")\n";
  }

  // Get the particle data
  auto ptData = data->GetPointData();
  auto ids = ptData->GetScalars("ID");
  auto types = ptData->GetScalars("Type");
  auto rads = ptData->GetVectors("Radius");
  auto axas = ptData->GetVectors("Axis a");
  auto axbs = ptData->GetVectors("Axis b");
  auto axcs = ptData->GetVectors("Axis c");
  auto vels = ptData->GetVectors("Velocity");
  auto omgs = ptData->GetVectors("Omega");
  auto fors = ptData->GetVectors("Force");
  auto moms = ptData->GetVectors("Moment");
  double scalar;
  double vector[3];
  for (int ii = 0; ii < numPts; ++ii) {
    ids->GetTuple(ii, &scalar);
    std::size_t id = static_cast<std::size_t>(scalar);

    types->GetTuple(ii, &scalar);
    int type = static_cast<int>(scalar);

    rads->GetTuple(ii, vector);
    Vec radius(vector[0], vector[1], vector[2]);

    axas->GetTuple(ii, vector);
    Vec axis_a(vector[0], vector[1], vector[2]);

    axbs->GetTuple(ii, vector);
    Vec axis_b(vector[0], vector[1], vector[2]);

    axcs->GetTuple(ii, vector);
    Vec axis_c(vector[0], vector[1], vector[2]);

    vels->GetTuple(ii, vector);
    Vec velocity(vector[0], vector[1], vector[2]);

    omgs->GetTuple(ii, vector);
    Vec omega(vector[0], vector[1], vector[2]);

    fors->GetTuple(ii, vector);
    Vec force(vector[0], vector[1], vector[2]);

    moms->GetTuple(ii, vector);
    Vec moment(vector[0], vector[1], vector[2]);

    //std::cout << id << ", ";
    //std::cout << type << ",";
    //std::cout << radius << ",";
    //std::cout << axis_a << "\n";
    //std::cout << force << "\n";

    DEMParticleP particle = std::make_shared<DEMParticle>(
      id, DEMParticle::DEMParticleShape::ELLIPSOID, 
      static_cast<DEMParticle::DEMParticleType>(type),
      radius, positions[ii], axis_a, axis_b, axis_c, 
      youngModulus, poissonRatio);

    particle->setPreviousVelocity(velocity);
    particle->setCurrentVelocity(velocity);
    particle->setPreviousAngularVelocity(omega);
    particle->setCurrentAngularVelocity(omega);
    particle->setForce(force);
    particle->setMoment(moment);

    particles.push_back(particle);
  }
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
    Ellip3D::IOUtil::decodeAndUncompress<T>(particleDataStr, numComp, output);
  if (!success) {
    std::cerr << "**ERROR** Could not decode and uncompress particle " << name
              << "\n";
    return false;
  }
  return true;
}

} // end namespace dem