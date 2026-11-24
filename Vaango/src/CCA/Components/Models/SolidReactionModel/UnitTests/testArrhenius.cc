#include <CCA/Components/Models/SolidReactionModel/Arrhenius.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Malloc/Allocator.h>

#include <gtest/gtest.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <cmath>

using namespace Uintah;

TEST(ArrheniusTest, ValueTest) {
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr root = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
  xmlDocSetRootElement(doc, root);
  
  xmlNodePtr node = xmlNewChild(root, nullptr, BAD_CAST "RateConstant", nullptr);
  xmlNewChild(node, nullptr, BAD_CAST "Ea", BAD_CAST "8314.462175");
  xmlNewChild(node, nullptr, BAD_CAST "A", BAD_CAST "1.0");

  ProblemSpecP ps = scinew ProblemSpec(node, false);
  
  Arrhenius model(ps);
  
  double R = 8.314462175;
  double T = 1000.0;
  double expected = 1.0 * std::exp(-8314.462175 / (R * T));
  
  EXPECT_NEAR(model.getConstant(T), expected, 1e-12);
  
  xmlFreeDoc(doc);
}

TEST(ArrheniusTest, ZeroEaTest) {
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr root = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
  xmlDocSetRootElement(doc, root);
  
  xmlNodePtr node = xmlNewChild(root, nullptr, BAD_CAST "RateConstant", nullptr);
  xmlNewChild(node, nullptr, BAD_CAST "Ea", BAD_CAST "0.0");
  xmlNewChild(node, nullptr, BAD_CAST "A", BAD_CAST "2.5");

  ProblemSpecP ps = scinew ProblemSpec(node, false);
  
  Arrhenius model(ps);
  
  EXPECT_DOUBLE_EQ(model.getConstant(300.0), 2.5);
  EXPECT_DOUBLE_EQ(model.getConstant(1000.0), 2.5);
  
  xmlFreeDoc(doc);
}
