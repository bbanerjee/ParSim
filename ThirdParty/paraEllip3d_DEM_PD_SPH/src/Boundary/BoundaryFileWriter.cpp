#include <Boundary/BoundaryFileWriter.h>
#include <Boundary/BoundaryContainers.h>
#include <Core/Geometry/Box.h>
#include <Core/Math/IntVec.h>
#include <Core/Const/const.h>

using namespace dem;

/**
 * Create the boundary file in the original CSV format
 */
void
BoundaryFileWriter::writeCSV(std::size_t boundaryNum,
                             const std::string& outputFileName, 
                             const Box& allContainer) const
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

  // Get the container limits and center
  Vec domainMin = allContainer.getMinCorner();
  Vec domainMax = allContainer.getMaxCorner();
  Vec domainCen = allContainer.getCenter();

  // Write the boundaryNum flag
  ofs << std::setw(OWID) << domainMin.x() 
      << std::setw(OWID) << domainMin.y() 
      << std::setw(OWID) << domainMin.z()
      << std::setw(OWID) << domainMax.x() 
      << std::setw(OWID) << domainMax.y() 
      << std::setw(OWID) << domainMax.z() << std::endl;
  ofs << std::endl;
  ofs << std::setw(OWID) << boundaryNum << std::endl;
  ofs << std::endl;

  // Get the boundary flag type
  BoundaryFlag flag = getBoundaryFlag(boundaryNum);

  // Various cases
  switch(flag) {

  case ONLY_BOTTOM_BOUNDARY:

    writeHeader(0u, ofs);
    writeBoundaryID(BoundaryID::ZMINUS, ofs);
    writeZMinusCSV(domainMin, domainMax, domainCen, ofs);
    writeEmptyLine(ofs);
    break;

  case NO_TOP_BOUNDARY:

    writeHeader(1u, ofs);
    writeBoundaryID(BoundaryID::XMINUS, ofs);
    writeXMinusCSV(domainMin, domainMax, domainCen, ofs);
    writeExtraZPlusCSV(domainMin, domainMax, domainCen, ofs);
    writeEmptyLine(ofs);

    writeHeader(1u, ofs);
    writeBoundaryID(BoundaryID::XPLUS, ofs);
    writeXPlusCSV(domainMin, domainMax, domainCen, ofs);
    writeExtraZPlusCSV(domainMin, domainMax, domainCen, ofs);
    writeEmptyLine(ofs);

    writeHeader(1u, ofs);
    writeBoundaryID(BoundaryID::YMINUS, ofs);
    writeYMinusCSV(domainMin, domainMax, domainCen, ofs);
    writeExtraZPlusCSV(domainMin, domainMax, domainCen, ofs);
    writeEmptyLine(ofs);

    writeHeader(1u, ofs);
    writeBoundaryID(BoundaryID::YPLUS, ofs);
    writeYPlusCSV(domainMin, domainMax, domainCen, ofs);
    writeExtraZPlusCSV(domainMin, domainMax, domainCen, ofs);
    writeEmptyLine(ofs);

    writeHeader(1u, ofs);
    writeBoundaryID(BoundaryID::ZMINUS, ofs);
    writeZMinusCSV(domainMin, domainMax, domainCen, ofs);
    writeExtraZPlusCSV(domainMin, domainMax, domainCen, ofs);
    writeEmptyLine(ofs);
    break;

  default:

    writeHeader(0u, ofs);
    writeBoundaryID(BoundaryID::XMINUS, ofs);
    writeXMinusCSV(domainMin, domainMax, domainCen, ofs);
    writeEmptyLine(ofs);

    writeHeader(0u, ofs);
    writeBoundaryID(BoundaryID::XPLUS, ofs);
    writeXPlusCSV(domainMin, domainMax, domainCen, ofs);
    writeEmptyLine(ofs);

    writeHeader(0u, ofs);
    writeBoundaryID(BoundaryID::YMINUS, ofs);
    writeYMinusCSV(domainMin, domainMax, domainCen, ofs);
    writeEmptyLine(ofs);

    writeHeader(0u, ofs);
    writeBoundaryID(BoundaryID::YPLUS, ofs);
    writeYPlusCSV(domainMin, domainMax, domainCen, ofs);
    writeEmptyLine(ofs);

    writeHeader(0u, ofs);
    writeBoundaryID(BoundaryID::ZMINUS, ofs);
    writeZMinusCSV(domainMin, domainMax, domainCen, ofs);
    writeEmptyLine(ofs);

    writeHeader(0u, ofs);
    writeBoundaryID(BoundaryID::ZPLUS, ofs);
    writeZPlusCSV(domainMin, domainMax, domainCen, ofs);
    writeEmptyLine(ofs);
    break;
  }

  ofs.close();
}

/**
 * Create a boundary file in XML format 
 */
void
BoundaryFileWriter::writeXML(std::size_t boundaryNum,
                             const std::string& outputFileName, 
                             const Box& allContainer) const
{
  // Get the container limits and center
  Vec domainMin = allContainer.getMinCorner();
  Vec domainMax = allContainer.getMaxCorner();
  Vec domainCen = allContainer.getCenter();

  // Create empty document
  zen::XmlDoc doc("Ellip3D_input");

  // Create a proxy output
  zen::XmlOut xml(doc);

  // Write the title
  std::string title = "Ellip3D boundary XML file";
  xml["Meta"]["title"](title);

  // Write the container dimensions
  std::ostringstream stream;
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.precision(dem::OPREC);
  stream << "[" << allContainer.getMinCorner().x() << ", "
                << allContainer.getMinCorner().y() << ", "
                << allContainer.getMinCorner().z() << "]";
  xml["Boundary"]["containerMin"](stream.str());

  stream.str("");
  stream << "[" << allContainer.getMaxCorner().x() << ", "
                << allContainer.getMaxCorner().y() << ", "
                << allContainer.getMaxCorner().z() << "]";
  xml["Boundary"]["containerMax"](stream.str());

  // Get the boundary flag type
  BoundaryFlag flag = getBoundaryFlag(boundaryNum);

  // Various cases
  switch(flag) {

  case ONLY_BOTTOM_BOUNDARY:

    {
    zen::XmlElement& element = xml["Boundary"].ref();
    zen::XmlElement& child = element.addChild("boundary");
    zen::XmlOut xml_child(child);
    writeHeader(0u, xml_child);
    writeBoundaryID(BoundaryID::ZMINUS, xml_child);
    writeZMinusXML(domainMin, domainMax, domainCen, xml_child);
    }
    break;

  case NO_TOP_BOUNDARY:

    {
    zen::XmlElement& element = xml["Boundary"].ref();
    zen::XmlElement& child = element.addChild("boundary");
    zen::XmlOut xml_child(child);
    writeHeader(1u, xml_child);
    writeBoundaryID(BoundaryID::XMINUS, xml_child);
    writeXMinusXML(domainMin, domainMax, domainCen, xml_child);
    writeExtraZPlusXML(domainMin, domainMax, domainCen, xml_child);
    }

    {
    zen::XmlElement& element = xml["Boundary"].ref();
    zen::XmlElement& child = element.addChild("boundary");
    zen::XmlOut xml_child(child);
    writeHeader(1u, xml_child);
    writeBoundaryID(BoundaryID::XPLUS, xml_child);
    writeXPlusXML(domainMin, domainMax, domainCen, xml_child);
    writeExtraZPlusXML(domainMin, domainMax, domainCen, xml_child);
    }

    {
    zen::XmlElement& element = xml["Boundary"].ref();
    zen::XmlElement& child = element.addChild("boundary");
    zen::XmlOut xml_child(child);
    writeHeader(1u, xml_child);
    writeBoundaryID(BoundaryID::YMINUS, xml_child);
    writeYMinusXML(domainMin, domainMax, domainCen, xml_child);
    writeExtraZPlusXML(domainMin, domainMax, domainCen, xml_child);
    }

    {
    zen::XmlElement& element = xml["Boundary"].ref();
    zen::XmlElement& child = element.addChild("boundary");
    zen::XmlOut xml_child(child);
    writeHeader(1u, xml_child);
    writeBoundaryID(BoundaryID::YPLUS, xml_child);
    writeYPlusXML(domainMin, domainMax, domainCen, xml_child);
    writeExtraZPlusXML(domainMin, domainMax, domainCen, xml_child);
    }

    {
    zen::XmlElement& element = xml["Boundary"].ref();
    zen::XmlElement& child = element.addChild("boundary");
    zen::XmlOut xml_child(child);
    writeHeader(1u, xml_child);
    writeBoundaryID(BoundaryID::ZMINUS, xml_child);
    writeZMinusXML(domainMin, domainMax, domainCen, xml_child);
    writeExtraZPlusXML(domainMin, domainMax, domainCen, xml_child);
    }
    break;

  default:

    {
    zen::XmlElement& element = xml["Boundary"].ref();
    zen::XmlElement& child = element.addChild("boundary");
    zen::XmlOut xml_child(child);
    writeHeader(0u, xml_child);
    writeBoundaryID(BoundaryID::XMINUS, xml_child);
    writeXMinusXML(domainMin, domainMax, domainCen, xml_child);
    }

    {
    zen::XmlElement& element = xml["Boundary"].ref();
    zen::XmlElement& child = element.addChild("boundary");
    zen::XmlOut xml_child(child);
    writeHeader(0u, xml_child);
    writeBoundaryID(BoundaryID::XPLUS, xml_child);
    writeXPlusXML(domainMin, domainMax, domainCen, xml_child);
    }

    {
    zen::XmlElement& element = xml["Boundary"].ref();
    zen::XmlElement& child = element.addChild("boundary");
    zen::XmlOut xml_child(child);
    writeHeader(0u, xml_child);
    writeBoundaryID(BoundaryID::YMINUS, xml_child);
    writeYMinusXML(domainMin, domainMax, domainCen, xml_child);
    }

    {
    zen::XmlElement& element = xml["Boundary"].ref();
    zen::XmlElement& child = element.addChild("boundary");
    zen::XmlOut xml_child(child);
    writeHeader(0u, xml_child);
    writeBoundaryID(BoundaryID::YPLUS, xml_child);
    writeYPlusXML(domainMin, domainMax, domainCen, xml_child);
    }

    {
    zen::XmlElement& element = xml["Boundary"].ref();
    zen::XmlElement& child = element.addChild("boundary");
    zen::XmlOut xml_child(child);
    writeHeader(0u, xml_child);
    writeBoundaryID(BoundaryID::ZMINUS, xml_child);
    writeZMinusXML(domainMin, domainMax, domainCen, xml_child);
    }

    {
    zen::XmlElement& element = xml["Boundary"].ref();
    zen::XmlElement& child = element.addChild("boundary");
    zen::XmlOut xml_child(child);
    writeHeader(0u, xml_child);
    writeBoundaryID(BoundaryID::ZPLUS, xml_child);
    writeZPlusXML(domainMin, domainMax, domainCen, xml_child);
    }
    break;
  }

  try {
    zen::save(doc, outputFileName);
  } catch (const zen::XmlFileError& err) {
    std::cerr << "**ERROR**: Could not write XML boundary file "
              << outputFileName << " for writing : in " 
              << __FILE__ << ":" << __LINE__
              << std::endl;
    exit(-1);
  }
}

/**
 * For XML output
 */
void
BoundaryFileWriter::writeHeader(std::size_t extraBoundaries,
                                std::ofstream& ofs) const
{
  ofs << std::setw(OWID) << 1 
      << std::setw(OWID) << extraBoundaries << std::endl;
}

void
BoundaryFileWriter::writeBoundaryID(BoundaryID id,
                                    std::ofstream& ofs) const
{
  if (id == BoundaryID::NONE) {
    ofs << std::setw(OWID) << " ";
  } else {
    ofs << std::setw(OWID) << id;
  }
}

void
BoundaryFileWriter::writeEmptyLine(std::ofstream& ofs)  const
{
  ofs << std::endl;
}

void
BoundaryFileWriter::writeExtraZPlusCSV(const Vec& domainMin, const Vec& domainMax,
                                       const Vec& domainCen, std::ofstream& ofs) const
{
  writeBoundaryID(BoundaryID::NONE, ofs);
  writeZPlusCSV(domainMin, domainMax, domainCen, ofs);
}

void
BoundaryFileWriter::writeXMinusCSV(const Vec& domainMin, const Vec& domainMax,
                                   const Vec& domainCen, std::ofstream& ofs) const
{
  ofs << std::setw(OWID) << -1 << std::setw(OWID) << 0 << std::setw(OWID) << 0
      << std::setw(OWID) << domainMin.x() 
      << std::setw(OWID) << domainCen.y() 
      << std::setw(OWID) << domainCen.z() << std::endl;
}

void
BoundaryFileWriter::writeXPlusCSV(const Vec& domainMin, const Vec& domainMax,
                                  const Vec& domainCen, std::ofstream& ofs) const
{
  ofs << std::setw(OWID) << 1 << std::setw(OWID) << 0 << std::setw(OWID) << 0
      << std::setw(OWID) << domainMax.x() 
      << std::setw(OWID) << domainCen.y() 
      << std::setw(OWID) << domainCen.z() << std::endl;
}

void
BoundaryFileWriter::writeYMinusCSV(const Vec& domainMin, const Vec& domainMax,
                                   const Vec& domainCen, std::ofstream& ofs) const
{
  ofs << std::setw(OWID) << 0 << std::setw(OWID) << -1 << std::setw(OWID) << 0
      << std::setw(OWID) << domainCen.x() 
      << std::setw(OWID) << domainMin.y() 
      << std::setw(OWID) << domainCen.z() << std::endl;
}

void
BoundaryFileWriter::writeYPlusCSV(const Vec& domainMin, const Vec& domainMax,
                                  const Vec& domainCen, std::ofstream& ofs) const
{
  ofs << std::setw(OWID) << 0 << std::setw(OWID) << 1 << std::setw(OWID) << 0
      << std::setw(OWID) << domainCen.x() 
      << std::setw(OWID) << domainMax.y() 
      << std::setw(OWID) << domainCen.z() << std::endl;
}

void
BoundaryFileWriter::writeZMinusCSV(const Vec& domainMin, const Vec& domainMax,
                                   const Vec& domainCen, std::ofstream& ofs) const
{
  ofs << std::setw(OWID) << 0 << std::setw(OWID) << 0 << std::setw(OWID) << -1
      << std::setw(OWID) << domainCen.x() 
      << std::setw(OWID) << domainCen.y() 
      << std::setw(OWID) << domainMin.z() << std::endl;
}

void
BoundaryFileWriter::writeZPlusCSV(const Vec& domainMin, const Vec& domainMax,
                                  const Vec& domainCen, std::ofstream& ofs) const
{
  ofs << std::setw(OWID) << 0 << std::setw(OWID) << 0 << std::setw(OWID) << 1
      << std::setw(OWID) << domainCen.x() 
      << std::setw(OWID) << domainCen.y() 
      << std::setw(OWID) << domainMax.z() << std::endl;
}

/**
 * For XML output
 */
void
BoundaryFileWriter::writeHeader(std::size_t extraBoundaries,
                                zen::XmlOut& xml) const
{
  xml.attribute("type", "plane");
}

void
BoundaryFileWriter::writeBoundaryID(BoundaryID id,
                                    zen::XmlOut& xml) const
{
  if (id != BoundaryID::NONE) {
    xml.attribute("id", static_cast<int>(id));
  }
}

void
BoundaryFileWriter::writeExtraZPlusXML(const Vec& domainMin, const Vec& domainMax,
                                       const Vec& domainCen, zen::XmlOut& xml) const
{
  zen::XmlElement& element = xml.ref();
  zen::XmlElement& child = element.addChild("extraEdge");
  zen::XmlOut xml_extra(child);
  writeZPlusXML(domainMin, domainMax, domainCen, xml_extra);
}

void
BoundaryFileWriter::writeXMinusXML(const Vec& domainMin, const Vec& domainMax,
                                   const Vec& domainCen, zen::XmlOut& xml) const
{
  xml["direction"]("[-1, 0, 0]");
  std::ostringstream stream;
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.precision(dem::OPREC);
  stream << "[" 
         << domainMin.x() << ", " << domainCen.y() << ", " << domainCen.z() 
         << "]";
  xml["position"](stream.str());
}

void
BoundaryFileWriter::writeXPlusXML(const Vec& domainMin, const Vec& domainMax,
                                  const Vec& domainCen, zen::XmlOut& xml) const
{
  xml["direction"]("[1, 0, 0]");
  std::ostringstream stream;
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.precision(dem::OPREC);
  stream << "[" 
         << domainMax.x() << ", " << domainCen.y() << ", " << domainCen.z() 
         << "]";
  xml["position"](stream.str());
}

void
BoundaryFileWriter::writeYMinusXML(const Vec& domainMin, const Vec& domainMax,
                                   const Vec& domainCen, zen::XmlOut& xml) const
{
  xml["direction"]("[0, -1, 0]");
  std::ostringstream stream;
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.precision(dem::OPREC);
  stream << "[" 
         << domainCen.x() << ", " << domainMin.y() << ", " << domainCen.z()
         << "]";
  xml["position"](stream.str());
}

void
BoundaryFileWriter::writeYPlusXML(const Vec& domainMin, const Vec& domainMax,
                                  const Vec& domainCen, zen::XmlOut& xml) const
{
  xml["direction"]("[0, 1, 0]");
  std::ostringstream stream;
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.precision(dem::OPREC);
  stream << "[" 
         << domainCen.x() << ", " << domainMax.y() << ", " << domainCen.z()
         << "]";
  xml["position"](stream.str());
}

void
BoundaryFileWriter::writeZMinusXML(const Vec& domainMin, const Vec& domainMax,
                                   const Vec& domainCen, zen::XmlOut& xml) const
{
  xml["direction"]("[0, 0, -1]");
  std::ostringstream stream;
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.precision(dem::OPREC);
  stream << "[" 
         << domainCen.x() << ", " << domainCen.y() << ", " << domainMin.z()
         << "]";
  xml["position"](stream.str());
}

void
BoundaryFileWriter::writeZPlusXML(const Vec& domainMin, const Vec& domainMax,
                                  const Vec& domainCen, zen::XmlOut& xml) const
{
  xml["direction"]("[0, 0, 1]");
  std::ostringstream stream;
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.precision(dem::OPREC);
  stream << "[" 
         << domainCen.x() << ", " << domainCen.y() << ", " << domainMax.z()
         << "]";
  xml["position"](stream.str());
}
