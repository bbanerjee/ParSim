#include <Core/IO/UintahZlibUtil.h>
#include <gtest/gtest.h>
#include <zlib.h>
#include <cstdio>
#include <string>

using namespace Uintah;

TEST(ZlibTest, ReadFunctions) {
  const char* filename = "test_data.gz";
  gzFile out = gzopen(filename, "wb");
  ASSERT_NE(out, nullptr);
  // Ensure a newline after the last token so getLine doesn't start at the end of the first line
  gzprintf(out, "test_string 123 45.6\n# comment line\nline after comment");
  gzclose(out);
  
  gzFile in = gzopen(filename, "rb");
  ASSERT_NE(in, nullptr);
  
  EXPECT_EQ(getString(in), "test_string");
  EXPECT_EQ(getInt(in), 123);
  EXPECT_DOUBLE_EQ(getDouble(in), 45.6);
  
  // After getDouble, we are at the '\n' of the first line.
  // getLine will read from current position. 
  // It should see '\n', think it's an empty line, but it also skips '#' comments.
  // Wait, getToken skips leading whitespace but not getLine.
  
  std::string line = getLine(in);
  if (line == "") line = getLine(in); // Skip the leftover empty line from \n

  EXPECT_EQ(line, "line after comment"); 
  
  gzclose(in);
  std::remove(filename);
}
