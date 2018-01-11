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

#include <InputOutput/DEMParticleContactFileReaderXML.h>
#include <DiscreteElements/DEMParticle.h>
#include <DiscreteElements/DEMContact.h>
#include <Core/Const/Constants.h>
#include <Core/Types/IntegerTypes.h>
#include <iomanip>
#include <exception>

using namespace dem;

DEMParticleContactFileReaderXML::DEMParticleContactFileReaderXML(const std::string& inputFileName) 
{
  // Read the input file
  try {
    //std::cout << "Input file name= " << inputFileName << "\n";
    d_doc = zen::load(inputFileName);
  } catch (const zen::XmlFileError& err) {
    std::cerr << "*ERROR** Could not read input file " << inputFileName << "\n";
    std::cerr << "    Error # = " << err.lastError << "\n";
    throw false;
  } catch (const zen::XmlParsingError& err) {
    std::cerr << "*ERROR** Could not read input file " << inputFileName << "\n";
    std::cerr << "    Parse Error in line: " << err.row + 1
              << " col: " << err.col << "\n";
    throw false;
  }

  // Check whether this is the right type of input file
  if (d_doc.root().getNameAs<std::string>() != "Ellip3D_input") {
    std::cerr << "*ERROR** Could not find tag <Ellip3D_input> in input file "
              << inputFileName << "\n";
    throw false;
  }
}

DEMParticleContactFileReaderXML::~DEMParticleContactFileReaderXML() 
{
}

/**
 * Read the particle contacts in XML format
 */
void
DEMParticleContactFileReaderXML::read(const DEMParticlePArray& particles,
                                      DEMContactArray& contacts)
{
  // Create a binary tree with particle IDs as key and particle pointer as value
  std::map<ParticleID, DEMParticle*> particleRP;
  for (const auto& particle : particles) {
    particleRP[particle->getId()] = particle.get();
  }
  
  // Load the document into input proxy for easier element access
  zen::XmlIn xml(d_doc);

  // Loop thru processes in the input file
  for (auto process_ps = xml["Process"]; process_ps; process_ps.next()) {

    // Loop thru contacts in each process
    for (auto contact_ps = process_ps["contact"]; contact_ps; contact_ps.next()) {

      // Get the attributes
      ParticleID particleID1 = 0;
      ParticleID particleID2 = 0;
      contact_ps.attribute("particle1", particleID1);
      contact_ps.attribute("particle2", particleID2);

      // Read the contact data
      DEMContactData data;
      data.read(contact_ps);
      
      // Create a DEMContact object
      contacts.push_back(DEMContact(particleRP[particleID1], particleRP[particleID2],
                                    data));
    }
  }
}
