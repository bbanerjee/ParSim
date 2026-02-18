#include <CCA/Components/Models/SolidReactionModel/NthOrderModel.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Malloc/Allocator.h>

#include <gtest/gtest.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <cmath>

using namespace Uintah;

TEST(NthOrderModelTest, ValueTest) {
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr root = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
  xmlDocSetRootElement(doc, root);
  
  xmlNodePtr node = xmlNewChild(root, nullptr, BAD_CAST "RateModel", nullptr);
  xmlNewChild(node, nullptr, BAD_CAST "n", BAD_CAST "2");

  ProblemSpecP ps = scinew ProblemSpec(node, false);
  
  NthOrderModel model(ps);
  
  double f = 0.3;
  double expected = std::pow(1.0 - f, 2.0);
  
  EXPECT_NEAR(model.getDifferentialFractionChange(f), expected, 1e-12);
  
  xmlFreeDoc(doc);
}

TEST(NthOrderModelTest, ZeroOrderTest) {
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr root = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
  xmlDocSetRootElement(doc, root);
  
  xmlNodePtr node = xmlNewChild(root, nullptr, BAD_CAST "RateModel", nullptr);
  xmlNewChild(node, nullptr, BAD_CAST "n", BAD_CAST "0");

  ProblemSpecP ps = scinew ProblemSpec(node, false);
  
  NthOrderModel model(ps);
  
  EXPECT_DOUBLE_EQ(model.getDifferentialFractionChange(0.5), 1.0);
  
  xmlFreeDoc(doc);
}
