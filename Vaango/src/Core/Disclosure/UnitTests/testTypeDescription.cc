#include <Core/Disclosure/TypeDescription.h>
#include <Core/Disclosure/TypeUtils.h>
#include <gtest/gtest.h>
#include <string>

using namespace Uintah;

// Note: TypeDescription uses a global map that is populated when 
// fun_getTypeDescription is called for various types.

TEST(TypeDescriptionTest, BasicTypes) {
  // We must ensure the types are registered by calling the fun_ versions
  double d;
  fun_getTypeDescription(&d);
  
  const TypeDescription* td = TypeDescription::lookupType("double");
  ASSERT_NE(td, nullptr);
  EXPECT_EQ(td->getName(), "double");
  EXPECT_EQ(td->getType(), TypeDescription::Type::double_type);
  EXPECT_TRUE(td->isFlat());
}

TEST(TypeDescriptionTest, IntType) {
  int i;
  fun_getTypeDescription(&i);

  const TypeDescription* td = TypeDescription::lookupType("int");
  ASSERT_NE(td, nullptr);
  EXPECT_EQ(td->getName(), "int");
  EXPECT_EQ(td->getType(), TypeDescription::Type::int_type);
}

TEST(TypeDescriptionTest, VectorType) {
  Vector v;
  fun_getTypeDescription(&v);

  const TypeDescription* td = TypeDescription::lookupType("Vector");
  ASSERT_NE(td, nullptr);
  EXPECT_EQ(td->getName(), "Vector");
  EXPECT_EQ(td->getType(), TypeDescription::Type::Vector);
}

TEST(TypeDescriptionTest, ToString) {
  EXPECT_EQ(TypeDescription::toString(TypeDescription::Type::double_type), "double_type");
  EXPECT_EQ(TypeDescription::toString(TypeDescription::Type::int_type), "int_type");
  EXPECT_EQ(TypeDescription::toString(TypeDescription::Type::Vector), "Vector");
}

TEST(TypeUtilsTest, FunGetDescription) {
  double d;
  const TypeDescription* td = fun_getTypeDescription(&d);
  ASSERT_NE(td, nullptr);
  EXPECT_EQ(td->getName(), "double");
  
  int i;
  td = fun_getTypeDescription(&i);
  ASSERT_NE(td, nullptr);
  EXPECT_EQ(td->getName(), "int");
}
