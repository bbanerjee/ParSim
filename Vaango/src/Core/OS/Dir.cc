/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */
#include <Core/OS/Dir.h>

#include <Core/Exceptions/ErrnoException.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Util/FileUtils.h>

#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

namespace Uintah {

Dir
Dir::create(const std::string& name)
{
  try {
    fs::create_directories(name);
    fs::permissions(name,
                    fs::perms::owner_all | fs::perms::group_all,
                    fs::perm_options::add);
  } catch (fs::filesystem_error& err) {
    throw InternalError("Dir::could not create directory: " + name + " " +
                          err.what(),
                        __FILE__,
                        __LINE__);
  }
  return Dir(name);
}

Dir::Dir(const std::string& name)
  : d_name(name)
{
}

Dir::Dir(const Dir& dir)
  : d_name(dir.d_name)
{
}

Dir&
Dir::operator=(const Dir& copy)
{
  d_name = copy.d_name;
  return *this;
}

void
Dir::remove(bool throwOnError)
{
  try {
    fs::remove(d_name);
  } catch (fs::filesystem_error& err) {
    InternalError exception("Dir::remove()::rmdir: " + d_name + " " +
                              err.what(),
                            __FILE__,
                            __LINE__);
    if (throwOnError) {
      throw exception;
    } else {
      std::cerr << "WARNING: " << exception.message() << "\n";
    }
  }
  return;
}

// Removes a directory (and all its contents).  Returns true if
// it succeeds.
bool
Dir::removeDir(const char* dirName)
{
  try {
    fs::remove_all(dirName);
  } catch (fs::filesystem_error& err) {
    std::cerr << "Error, rmdir failed for '" << dirName << "'\n";
    std::cerr << "  err is " << err.what() << "\n";
    return false;
  }
  return true;
}

void
Dir::forceRemove(bool throwOnError)
{
  try {
    fs::remove_all(d_name);
  } catch (fs::filesystem_error& err) {
    InternalError exception("Dir::remove()::rmdir: " + d_name + " " +
                              err.what(),
                            __FILE__,
                            __LINE__);
    if (throwOnError) {
      throw exception;
    } else {
      std::cerr << "WARNING: " << exception.message() << "\n";
    }
  }
}

void
Dir::remove(const std::string& filename, bool throwOnError)
{
  if (!exists()) {
    return;
  }
  std::string filepath = d_name + "/" + filename;
  try {
    fs::remove(filepath);
  } catch (fs::filesystem_error& err) {
    InternalError exception("Dir::remove()::rmdir: " + d_name + " " +
                              err.what(),
                            __FILE__,
                            __LINE__);
    if (throwOnError) {
      throw exception;
    } else {
      std::cerr << "WARNING: " << exception.message() << "\n";
    }
  }
  return;
}

Dir
Dir::createSubdir(const std::string& sub)
{
  return create(d_name + "/" + sub);
}

Dir
Dir::getSubdir(const std::string& sub)
{
  // This should probably do more
  return Dir(d_name + "/" + sub);
}

void
Dir::copy(Dir& destDir)
{
  const auto copyOptions = fs::copy_options::recursive;
  try {
    fs::copy(d_name, destDir.d_name, copyOptions);
  } catch (fs::filesystem_error& err) {
    throw InternalError(std::string("Dir::copy failed to copy: ") + d_name,
                        __FILE__,
                        __LINE__);
  }
  return;
}

void
Dir::move(Dir& destDir)
{
  try {
    fs::rename(d_name, destDir.d_name);
  } catch (fs::filesystem_error& err) {
    throw InternalError(std::string("Dir::move failed to move: ") + d_name,
                        __FILE__,
                        __LINE__);
  }
  return;
}

void
Dir::copy(const std::string& filename, Dir& destDir)
{
  std::string filepath = d_name + "/" + filename;
  try {
    fs::copy(filepath, destDir.d_name + "/" + filename);
  } catch (fs::filesystem_error& err) {
    throw InternalError(std::string("Dir::copy failed to copy: ") + filepath,
                        __FILE__,
                        __LINE__);
  }
  return;
}

void
Dir::move(const std::string& filename, Dir& destDir)
{
  std::string filepath = d_name + "/" + filename;
  try {
    fs::rename(filepath, destDir.d_name + "/" + filename);
  } catch (fs::filesystem_error& err) {
    throw InternalError(std::string("Dir::move failed to move: ") + filepath,
                        __FILE__,
                        __LINE__);
  }
  return;
}

void
Dir::getFilenamesBySuffix(const std::string& suffix,
                          std::vector<std::string>& filenames)
{
  if (!fs::exists(d_name)) {
    return;
  }

  std::string ext = "." + suffix;
  for (const auto& file : fs::directory_iterator{ d_name }) {
    if (file.path().extension() == ext) {
      filenames.push_back(file.path().string());
      std::cout << "  Found " << file.path().string() << "\n";
    }
  }
}

bool
Dir::exists()
{
  if (fs::exists(d_name)) {
    return true;
  }
  return false;
}

} // namespace Uintah