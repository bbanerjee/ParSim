#include <Core/DataArchive/DataArchive.h>
#include <Core/Exceptions/InternalError.h>
#include <gtest/gtest.h>
#include <string>
#include <fstream>
#include <filesystem>

using namespace Uintah;
namespace fs = std::filesystem;

TEST(DataArchiveTest, Construction) {
  std::string testUda = "test_data.uda";
  fs::create_directory(testUda);
  
  std::string indexXml = testUda + "/index.xml";
  std::ofstream ofs(indexXml);
  ofs << "<?xml version='1.0' encoding='ISO-8859-1' ?>\n";
  ofs << "<Uintah_update>\n";
  ofs << "  <outputFormat>UDA</outputFormat>\n";
  ofs << "</Uintah_update>\n";
  ofs.close();

  try {
    DataArchive da(testUda);
    EXPECT_EQ(da.name(), testUda);
    EXPECT_FALSE(da.isPIDXFormat());
  } catch (const InternalError& e) {
    FAIL() << "DataArchive construction failed with InternalError: " << e.message();
  } catch (...) {
    FAIL() << "DataArchive construction failed with unknown exception";
  }

  fs::remove_all(testUda);
}
