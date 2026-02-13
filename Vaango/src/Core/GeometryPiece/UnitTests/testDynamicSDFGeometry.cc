#include <Core/GeometryPiece/DynamicSDFGeometry.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Grid/Box.h>
#include <Core/Geometry/Point.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <gtest/gtest.h>

using namespace Uintah;

TEST(DynamicSDFGeometryTest, BasicCreation) {
  // Create a minimal XML for DynamicSDFGeometry
  // Note: TriGeometryPiece will try to read a file if we don't mock it well.
  // But for a unit test, we might need a dummy .dat file or a more complex mock.
  // For now, let's just test if we can at least check for null piece or catch expected exceptions.
  
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "dynamic_sdf");
  xmlDocSetRootElement(doc, rootNode);
  
  xmlNewChild(rootNode, nullptr, BAD_CAST "res", BAD_CAST "[10,10,10]");
  xmlNewChild(rootNode, nullptr, BAD_CAST "tri", BAD_CAST "");
  // We'll need a dummy file for TriGeometryPiece to not crash or throw if it actually tries to read.
  // Actually, TriGeometryPiece throws if file is missing.
  
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  
  // Since we don't have a real .dat file here, this constructor might throw.
  // We can at least check if it throws the right exception if we don't provide a file.
  EXPECT_THROW(scinew DynamicSDFGeometry(ps), Uintah::ProblemSetupException);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
