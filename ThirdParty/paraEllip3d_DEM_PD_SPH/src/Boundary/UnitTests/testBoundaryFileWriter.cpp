#include <Boundary/BoundaryFileWriter.h>
#include <Core/Util/Utility.h>
#include <gtest/gtest.h>

using namespace dem;

// Spherical particle - axis aligned
TEST(BoundaryFileWriterTest, writeCSV) {

  BoundaryFileWriter writer;

  std::size_t boundaryNum = 1u;
  Box spatialDomain(0, 0, 0, 1, 1, 1);

  std::string outputTxt = "test1.txt";
  writer.writeCSV(boundaryNum, outputTxt, spatialDomain);
  std::string outputXML= "test1.xml";
  writer.writeXML(boundaryNum, outputXML, spatialDomain);

  boundaryNum = 5u;
  outputTxt = "test2.txt";
  writer.writeCSV(boundaryNum, outputTxt, spatialDomain);
  outputXML= "test2.xml";
  writer.writeXML(boundaryNum, outputXML, spatialDomain);

  boundaryNum = 6u;
  outputTxt = "test3.txt";
  writer.writeCSV(boundaryNum, outputTxt, spatialDomain);
  outputXML= "test3.xml";
  writer.writeXML(boundaryNum, outputXML, spatialDomain);
}

