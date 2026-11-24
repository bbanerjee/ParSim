#include <CCA/Components/Models/SolidReactionModel/ModifiedArrhenius.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Malloc/Allocator.h>

#include <gtest/gtest.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <cmath>

using namespace Uintah;

TEST(ModifiedArrheniusTest, ValueTest) {
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr root = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
  xmlDocSetRootElement(doc, root);
  
  xmlNodePtr node = xmlNewChild(root, nullptr, BAD_CAST "RateConstant", nullptr);
  xmlNewChild(node, nullptr, BAD_CAST "Ea", BAD_CAST "8.314462175");
  xmlNewChild(node, nullptr, BAD_CAST "A", BAD_CAST "1.0");
  xmlNewChild(node, nullptr, BAD_CAST "b", BAD_CAST "2.0");

  ProblemSpecP ps = scinew ProblemSpec(node, false);
  
  ModifiedArrhenius model(ps);
  
  double R = 8.314462175;
  double T = 1.0;
  double expected = 1.0 * std::pow(T, 2.0) * std::exp(-8.314462175 / (R * T));
  
  EXPECT_NEAR(model.getConstant(T), expected, 1e-12);
  
  xmlFreeDoc(doc);
}

TEST(ModifiedArrheniusTest, ZeroEaTest) {
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr root = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
  xmlDocSetRootElement(doc, root);
  
  xmlNodePtr node = xmlNewChild(root, nullptr, BAD_CAST "RateConstant", nullptr);
  xmlNewChild(node, nullptr, BAD_CAST "Ea", BAD_CAST "0.0");
  xmlNewChild(node, nullptr, BAD_CAST "A", BAD_CAST "2.0");
  xmlNewChild(node, nullptr, BAD_CAST "b", BAD_CAST "0.5");

  ProblemSpecP ps = scinew ProblemSpec(node, false);
  
  ModifiedArrhenius model(ps);
  
  EXPECT_DOUBLE_EQ(model.getConstant(4.0), 4.0);
  EXPECT_DOUBLE_EQ(model.getConstant(16.0), 8.0);
  
  xmlFreeDoc(doc);
}
