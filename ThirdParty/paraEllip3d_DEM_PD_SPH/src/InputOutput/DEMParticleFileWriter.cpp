/*
 * The MIT License
 *
 * Copyright (c) 2017-2018 Parresia Research Limited, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <InputOutput/DEMParticleFileWriter.h>
#include <DiscreteElements/DEMParticle.h>
#include <Core/Const/Constants.h>
#include <Core/Math/Vec.h>
#include <iomanip>

using namespace dem;

/**
 * Create the particle input file in the original CSV format
 */
void
DEMParticleFileWriter::writeCSV(const DEMParticlePArray& particles,
                                const Gradation& gradation,
                                const std::string& outputFileName) const 
{
  // Open the output file
  std::ofstream ofs(outputFileName);
  if (!ofs) {
    std::cerr << "**ERROR**: Could not open boundary file "
              << outputFileName << " for writing : in " 
              << __FILE__ << ":" << __LINE__
              << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(dem::OPREC);
  ofs.width(dem::OWID);

  ofs << std::setw(OWID) << particles.size() << std::endl;
  ofs << std::setw(OWID) << "id" << std::setw(OWID) << "type" << std::setw(OWID)
      << "radius_a" << std::setw(OWID) << "radius_b" << std::setw(OWID)
      << "radius_c" << std::setw(OWID) << "position_x" << std::setw(OWID)
      << "position_y" << std::setw(OWID) << "position_z" << std::setw(OWID)
      << "axle_a_x" << std::setw(OWID) << "axle_a_y" << std::setw(OWID)
      << "axle_a_z" << std::setw(OWID) << "axle_b_x" << std::setw(OWID)
      << "axle_b_y" << std::setw(OWID) << "axle_b_z" << std::setw(OWID)
      << "axle_c_x" << std::setw(OWID) << "axle_c_y" << std::setw(OWID)
      << "axle_c_z" << std::setw(OWID) << "velocity_x" << std::setw(OWID)
      << "velocity_y" << std::setw(OWID) << "velocity_z" << std::setw(OWID)
      << "omga_x" << std::setw(OWID) << "omga_y" << std::setw(OWID) << "omga_z"
      << std::setw(OWID) << "force_x" << std::setw(OWID) << "force_y"
      << std::setw(OWID) << "force_z" << std::setw(OWID) << "moment_x"
      << std::setw(OWID) << "moment_y" << std::setw(OWID) << "moment_z"
      << std::endl;

  Vec vObj;
  for (const auto& part : particles) {
    ofs << std::setw(OWID) << part->getId() << std::setw(OWID)
        << part->getType() << std::setw(OWID) << part->getA() << std::setw(OWID)
        << part->getB() << std::setw(OWID) << part->getC();

    vObj = part->currentPosition();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->getCurrDirecA();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->getCurrDirecB();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->getCurrDirecC();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->currentVelocity();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->currentOmega();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->getForce();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->getMoment();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z() << std::endl;
  }

  // Write the gradation data
  gradation.outputCSV(ofs);

  ofs.close();
}

/**
 * Create a particle input file in XML format 
 */
void
DEMParticleFileWriter::writeXML(const DEMParticlePArray& particles,
                                const Gradation& gradation,
                                const std::string& outputFileName) const 
{
  // Create empty document
  zen::XmlDoc doc("Ellip3D_input");

  // Create a proxy output
  zen::XmlOut xml(doc);

  // Find the number of unique particle types
  std::vector<std::size_t> uniqueParticleTypes;
  for (const auto& particle : particles) {
    uniqueParticleTypes.push_back(particle->getType());
  }
  std::sort(uniqueParticleTypes.begin(), uniqueParticleTypes.end());
  auto last = std::unique(uniqueParticleTypes.begin(), uniqueParticleTypes.end());
  uniqueParticleTypes.erase(last, uniqueParticleTypes.end());

  // Loop through the particle types
  for (const auto particleType : uniqueParticleTypes) {

    // Loop through the particles and store in vectors
    std::size_t numParticles = particles.size();
    std::vector<std::size_t> particleIDs;
    std::vector<Vec> particleRadii;
    std::vector<Vec> particleAxleA, particleAxleB, particleAxleC;
    std::vector<Vec> particlePos, particleVel, particleRot;
    std::vector<Vec> particleForce, particleMoment;
    particleIDs.reserve(numParticles);
    particleRadii.reserve(numParticles);
    particleAxleA.reserve(numParticles); 
    particleAxleB.reserve(numParticles);
    particleAxleC.reserve(numParticles);
    particlePos.reserve(numParticles);
    particleVel.reserve(numParticles); 
    particleRot.reserve(numParticles);
    particleForce.reserve(numParticles);
    particleMoment.reserve(numParticles);
    for (const auto& particle : particles) {
      if (particle->getType() == particleType) {
        particleIDs.push_back(particle->getId());
        particleRadii.push_back(Vec(particle->getA(), particle->getB(), 
                                    particle->getC()));
        particleAxleA.push_back(particle->getCurrDirecA());
        particleAxleB.push_back(particle->getCurrDirecB());
        particleAxleC.push_back(particle->getCurrDirecC());
        particlePos.push_back(particle->currentPosition());
        particleVel.push_back(particle->currentVelocity());
        particleRot.push_back(particle->currentOmega());
        particleForce.push_back(particle->getForce());
        particleMoment.push_back(particle->getMoment());
      }
    }

    std::string type = DEMParticle::getDEMParticleShape(particleType);
    zen::XmlElement& root = xml.ref();
    zen::XmlElement& element = root.addChild("Particles");
    zen::XmlOut particles(element);
    particles.attribute("number", particleIDs.size());
    particles.attribute("type", type);
    particles.attribute("compression", "gzip");
    particles.attribute("encoding", "base64");

    // Write particle IDs
    writeParticleValues<std::size_t, 1>(element, "id", "none", particleIDs);

    // Write particle radii
    writeParticleValues<Vec, 3>(element, "radii", "mm", particleRadii);

    // Write particle axle_a, axle_b, axle_c
    writeParticleValues<Vec, 3>(element, "axle_a", "mm", particleAxleA);
    writeParticleValues<Vec, 3>(element, "axle_b", "mm", particleAxleB);
    writeParticleValues<Vec, 3>(element, "axle_c", "mm", particleAxleC);

    // Write particle position, velocity, omega
    writeParticleValues<Vec, 3>(element, "position", "m", particlePos);
    writeParticleValues<Vec, 3>(element, "velocity", "m/s", particleVel);
    writeParticleValues<Vec, 3>(element, "omega", "rad/s", particleRot);

    // Write particle force and moment
    writeParticleValues<Vec, 3>(element, "force", "N", particleForce);
    writeParticleValues<Vec, 3>(element, "moment", "N-m", particleMoment);
  }

  // Write the gradation data
  gradation.outputXML(xml);

  try {
    zen::save(doc, outputFileName);
  } catch (const zen::XmlFileError& err) {
    std::cerr << "**ERROR**: Could not write XML particle file "
              << outputFileName << " in " 
              << __FILE__ << ":" << __LINE__
              << std::endl;
    exit(-1);
  }
}

template <typename T, int numComponents>
void 
DEMParticleFileWriter::writeParticleValues(zen::XmlElement& element, 
                                           const std::string& name,
                                           const std::string& unit,
                                           std::vector<T>& particleData) const
{
  zen::XmlElement& child = element.addChild(name);
  zen::XmlOut xml_child(child);
  xml_child.attribute("unit", unit);
  xml_child.attribute("numComponents", numComponents);

  std::string str;
  if (!Ellip3D::Util::compressAndEncode<T>(particleData, numComponents, str)) {
    std::cerr << "**ERROR**: Could compresse and encode particle data for "
              << name << " in " 
              << __FILE__ << ":" << __LINE__
              << std::endl;
    exit(-1);
  }
  //std::ostringstream stream;
  //stream.setf(std::ios::scientific, std::ios::floatfield);
  //stream.precision(dem::OPREC);
  //for (const auto& value : particleData) {
  //  stream << value << " ";
  //}
  //xml_child(stream.str());
  xml_child(str);
}        