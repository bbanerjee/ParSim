#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Malloc/Allocator.h>
#include <gtest/gtest.h>
#include <string>

using namespace Uintah;

TEST(ProblemSpecTest, ParseBuffer) {
  std::string xml = R"(
    <Root>
      <Double>1.23</Double>
      <Int>456</Int>
      <String>hello</String>
      <Vector>[1, 2, 3]</Vector>
      <WithAttr attr="val">body</WithAttr>
    </Root>
  )";
    
  ProblemSpecP ps = scinew ProblemSpec(xml);
  ASSERT_TRUE(ps != nullptr);
  
  double d;
  ps->require("Double", d);
  EXPECT_DOUBLE_EQ(d, 1.23);
  
  int i;
  ps->require("Int", i);
  EXPECT_EQ(i, 456);
  
  std::string s;
  ps->require("String", s);
  EXPECT_EQ(s, "hello");
  
  Vector v;
  ps->require("Vector", v);
  EXPECT_EQ(v, Vector(1, 2, 3));
  
  ProblemSpecP child = ps->findBlock("WithAttr");
  ASSERT_TRUE(child != nullptr);
  std::string attr;
  child->getAttribute("attr", attr);
  EXPECT_EQ(attr, "val");
  EXPECT_EQ(child->getNodeValue(), "body");
}

TEST(ProblemSpecTest, GetWithDefault) {
  std::string xml = "<Root></Root>";
  ProblemSpecP ps = scinew ProblemSpec(xml);
  
  double d;
  ps->getWithDefault("Missing", d, 9.9);
  EXPECT_DOUBLE_EQ(d, 9.9);
}
