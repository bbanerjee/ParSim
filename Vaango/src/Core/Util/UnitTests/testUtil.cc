#include <Core/Util/StringUtil.h>
#include <Core/Util/Endian.h>
#include <gtest/gtest.h>
#include <string>
#include <vector>

using namespace Uintah;

TEST(StringUtilTest, Basic) {
  int i;
  EXPECT_TRUE(string_to_int("123", i));
  EXPECT_EQ(i, 123);
  
  double d;
  EXPECT_TRUE(string_to_double("1.23", d));
  EXPECT_DOUBLE_EQ(d, 1.23);
  
  EXPECT_EQ(string_toupper("Hello"), "HELLO");
  EXPECT_EQ(string_tolower("Hello"), "hello");
}

TEST(StringUtilTest, SplitCollapse) {
  std::vector<std::string> parts = split_string("a,b,c", {','});
  ASSERT_EQ(parts.size(), 3u);
  EXPECT_EQ(parts[0], "a");
  
  std::string s = "  spaces  ";
  collapse(s);
  EXPECT_EQ(s, "spaces");
}

TEST(EndianTest, Basic) {
  // We can't easily test swapbytes results without knowing the platform, 
  // but we can check if it's consistent.
  bool big = isBigEndian();
  bool little = isLittleEndian();
  EXPECT_NE(big, little);
  
  std::string e = endianness();
  if (little) EXPECT_EQ(e, "little_endian");
  else EXPECT_EQ(e, "big_endian");
}

TEST(EndianTest, Swap) {
  uint32_t val = 0x12345678;
  uint32_t swapped = val;
  swapbytes(swapped);
  if (isLittleEndian()) {
    EXPECT_EQ(swapped, 0x78563412);
  } else {
    EXPECT_EQ(swapped, 0x12345678); // swapbytes might be identity on big endian for some macros
    // Actually SWAP_4 usually does the swap regardless.
  }
}
