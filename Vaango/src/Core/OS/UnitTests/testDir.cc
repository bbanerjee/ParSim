#include <Core/Exceptions/InternalError.h>
#include <Core/OS/Dir.h>

#include <gtest/gtest.h>

#include <iostream>
#include <fstream>

using namespace Uintah;

TEST(DirTest, construction)
{

  Dir test_dir_init("test_dir");
  test_dir_init.forceRemove(true);

  Dir test_dir;
  Dir test_dir0("test_dir_0");
  Dir test_dir00(test_dir0);
  Dir test_dir01 = test_dir0;

  Dir test_dir1 = Dir::create("test_dir");
  Dir test_dir2 = Dir::create("test_dir");

  test_dir1.remove(true);
  test_dir2.remove(true);

  Dir test_dir3 = Dir::create("test_dir/first");
  Dir test_dir4 = Dir::create("test_dir/second");

  [[maybe_unused]] bool status = Dir::removeDir("./test_dir");

  Dir test_dir5 = Dir::create("test_dir/first");
  Dir test_dir6 = Dir::create("test_dir/second");
  Dir test_dir7("test_dir");
  test_dir7.forceRemove(true);

  Dir test_dir8 = Dir::create("test_dir/first");
  Dir test_dir9 = Dir::create("test_dir/second");
  test_dir7.remove("first", true);
  test_dir7.forceRemove(true);

  test_dir7.createSubdir("first");
  Dir test_dir10 = test_dir7.getSubdir("first");
  test_dir7.forceRemove(true);

  Dir test_dir11 = test_dir7.createSubdir("first");
  Dir test_dir12("test_dir/second");
  Dir test_dir13("test_dir/third");
  try {
    test_dir11.copy(test_dir12);
  } catch (const InternalError& err) {
    std::cout << "Copy error";
    throw;
  }
  try {
    test_dir11.move(test_dir13);
  } catch (const InternalError& err) {
    std::cout << "Move error";
    throw;
  }

  std::ofstream file1("test_dir/second/file1.txt");
  std::ofstream file2("test_dir/second/file2.txt");

  Dir test_dir16("test_dir/third");
  try {
    test_dir12.copy("file1.txt", test_dir16);
  } catch (const InternalError& err) {
    std::cout << "File copy error" << err.message();
    throw;
  }
  try {
    test_dir12.move("file2.txt", test_dir16);
  } catch (const InternalError& err) {
    std::cout << "File move error" << err.message();
    throw;
  }

  std::vector<std::string> filenames;
  test_dir16.getFilenamesBySuffix("txt", filenames);
  test_dir7.forceRemove(true);
}
