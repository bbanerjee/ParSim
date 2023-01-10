#include <Core/OS/Dir.h>

#include <gtest/gtest.h>

#include <iostream>

using namespace Uintah;

TEST(DirTest, construction) {

  Dir test_dir;
  Dir test_dir0("test_dir_0");
  Dir test_dir00(test_dir0);
  Dir test_dir01 = test_dir0;

  Dir test_dir1 = Dir::create("./test_dir");
  Dir test_dir2 = Dir::create("./test_dir");

  test_dir1.remove(true);
  test_dir2.remove(true);

  Dir test_dir3 = Dir::create("./test_dir/first");
  Dir test_dir4 = Dir::create("./test_dir/second");

  [[maybe_unused]] bool status = Dir::removeDir("./test_dir");

  Dir test_dir5 = Dir::create("./test_dir/first");
  Dir test_dir6 = Dir::create("./test_dir/second");
  Dir test_dir7("./test_dir");
  test_dir7.forceRemove(true);

  Dir test_dir8 = Dir::create("./test_dir/first");
  Dir test_dir9 = Dir::create("./test_dir/second");
  test_dir7.remove("first", true);
  test_dir7.forceRemove(true);

  test_dir7.createSubdir("first");
  Dir test_dir10 = test_dir7.getSubdir("first");
  test_dir7.forceRemove(true);

  Dir test_dir11 = test_dir7.createSubdir("first");
  Dir test_dir12("./test_dir/second");
  Dir test_dir13("./test_dir/third");
  test_dir11.copy(test_dir12);
  test_dir11.move(test_dir13);

  Dir test_dir14 = Dir::create("./test_dir/second/file1.txt");
  Dir test_dir15 = Dir::create("./test_dir/second/file2.txt");
  Dir test_dir16("./test_dir/third");
  test_dir12.copy("file1.txt", test_dir16);
  test_dir12.move("file2.txt", test_dir16);

  std::vector<std::string> filenames;
  test_dir16.getFilenamesBySuffix("txt", filenames);
  test_dir7.forceRemove(true);
}

