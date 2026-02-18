#include <CCA/Components/Models/SolidReactionModel/PowerModel.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Malloc/Allocator.h>

#include <gtest/gtest.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <cmath>

using namespace Uintah;

TEST(PowerModelTest, ValueTest) {
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr root = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
  xmlDocSetRootElement(doc, root);
  
  xmlNodePtr node = xmlNewChild(root, nullptr, BAD_CAST "RateModel", nullptr);
  xmlNewChild(node, nullptr, BAD_CAST "a", BAD_CAST "2.0");
  xmlNewChild(node, nullptr, BAD_CAST "b", BAD_CAST "0.5");

  ProblemSpecP ps = scinew ProblemSpec(node, false);
  
  PowerModel model(ps);
  
  EXPECT_DOUBLE_EQ(model.getDifferentialFractionChange(4.0), 4.0);
  EXPECT_DOUBLE_EQ(model.getDifferentialFractionChange(16.0), 8.0);
  
  xmlFreeDoc(doc);
}
