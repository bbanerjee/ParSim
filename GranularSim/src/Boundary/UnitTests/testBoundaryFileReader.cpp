#include <Boundary/BoundaryFileReader.h>
#include <Core/Geometry/OrientedBox.h>
#include <Core/Util/Utility.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(BoundaryFileReaderTest, readVTK) {

  BoundaryFileReader reader;

  std::string inputFile = "oriented_domain_00000.vtu";
  OrientedBox spatialDomain(Box(Vec(0,0,0), Vec(1,1,1)));
  reader.readVTK(inputFile, spatialDomain);
}

