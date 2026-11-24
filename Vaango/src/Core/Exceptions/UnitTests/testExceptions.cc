#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/FileNotFound.h>
#include <Core/Exceptions/ArrayIndexOutOfBounds.h>
#include <gtest/gtest.h>
#include <string>

using namespace Uintah;

TEST(ExceptionTest, InternalError) {
  try {
    throw InternalError("test error", __FILE__, __LINE__);
  } catch (const InternalError& e) {
    std::string msg = e.message();
    EXPECT_TRUE(msg.find("test error") != std::string::npos);
    EXPECT_STREQ(e.type(), "InternalError");
  }
}

TEST(ExceptionTest, ProblemSetupException) {
  try {
    throw ProblemSetupException("setup error", __FILE__, __LINE__);
  } catch (const ProblemSetupException& e) {
    std::string msg = e.message();
    EXPECT_TRUE(msg.find("setup error") != std::string::npos);
    EXPECT_STREQ(e.type(), "Uintah::Exceptions::ProblemSetupException");
  }
}

TEST(ExceptionTest, FileNotFound) {
  try {
    throw FileNotFound("missing_file.txt", __FILE__, __LINE__);
  } catch (const FileNotFound& e) {
    std::string msg = e.message();
    EXPECT_TRUE(msg.find("missing_file.txt") != std::string::npos);
    EXPECT_STREQ(e.type(), "FileNotFound");
  }
}

TEST(ExceptionTest, ArrayIndexOutOfBounds) {
  try {
    throw ArrayIndexOutOfBounds(10, 0, 5, __FILE__, __LINE__);
  } catch (const ArrayIndexOutOfBounds& e) {
    std::string msg = e.message();
    EXPECT_TRUE(msg.find("index 10") != std::string::npos);
    EXPECT_STREQ(e.type(), "ArrayIndexOutOfBounds");
  }
}
