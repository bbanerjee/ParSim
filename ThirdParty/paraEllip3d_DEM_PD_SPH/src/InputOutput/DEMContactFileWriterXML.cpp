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

#include <InputOutput/DEMContactFileWriterXML.h>
#include <DiscreteElements/DEMParticle.h>
#include <Core/Const/Constants.h>
#include <Core/Math/Vec.h>
#include <iomanip>

using namespace dem;

DEMContactFileWriterXML::DEMContactFileWriterXML(const std::string& outputFileName) 
  : d_doc("Ellip3D_input"),
    d_outputFileName(outputFileName+".xml")
{
}

DEMContactFileWriterXML::~DEMContactFileWriterXML() 
{
  try {
    zen::save(d_doc, d_outputFileName);
  } catch (const zen::XmlFileError& err) {
    std::cerr << "**ERROR**: Could not write XML boundary file "
              << d_outputFileName << " for writing : in "
              << __FILE__ << ":" << __LINE__
              << std::endl;
    exit(-1);
  }
}

/**
 * Write the boundary contacts in XML format
 */
void
DEMContactFileWriterXML::writeBoundaryContacts(const BoundaryPArray& boundaries)
{
  // Create a proxy output
  zen::XmlOut xml(d_doc);

  // Write the title
  std::string title = "Ellip3D boundary contact XML file";
  xml["Meta"]["title"](title);

  // Write the boundaries
  for (const auto& boundary : boundaries) {

    int boundaryID = static_cast<int>(boundary->getId());
    const BoundaryContactArray& contacts = boundary->getBoundaryContacts();

    zen::XmlElement& root = xml.ref();
    zen::XmlElement& child = root.addChild("boundary");
    zen::XmlOut xml_child(child);
    xml_child.attribute("id", boundaryID);
    xml_child.attribute("num_contacts", contacts.size());

    for (const auto& contact : contacts) {
      zen::XmlElement& child_contact = child.addChild("contact");
      zen::XmlOut xml_child_contact(child_contact);
      contact.write(xml_child_contact);
    }
  }
}

/**
 * Write the particle contacts in XML format
 */
void
DEMContactFileWriterXML::writeParticleContacts(const DEMParticlePArray& particles)
{

}