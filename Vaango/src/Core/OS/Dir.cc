/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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
#include <Core/Parallel/Parallel.h>
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
    fs::permissions(
      name, fs::perms::owner_all | fs::perms::group_all, fs::perm_options::add);
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
    InternalError exception(
      "Dir::remove()::rmdir: " + d_name + " " + err.what(), __FILE__, __LINE__);
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
    InternalError exception(
      "Dir::remove()::rmdir: " + d_name + " " + err.what(), __FILE__, __LINE__);
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
    InternalError exception(
      "Dir::remove()::rmdir: " + d_name + " " + err.what(), __FILE__, __LINE__);
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

// This version of createSubdir tries multiple times and
// throws an exception if it fails.
Dir
Dir::createSubdirPlus(const std::string& sub)
{
  Dir myDir;

  bool done = false;
  int tries = 0;

  while (!done) {

    try {
      tries++;

      if (tries > 500) {
        std::ostringstream warn;
        warn
          << " ERROR: Dir::createSubdirPlus() failed to create the directory ("
          << sub << ") after " << tries << " attempts.";
        throw InternalError(warn.str(), __FILE__, __LINE__);
      }

      myDir = createSubdir(sub);
      done  = true;
    } catch (ErrnoException& e) {
      if (e.getErrno() == EEXIST) {
        done  = true;
        myDir = getSubdir(sub);
      }
    }
  }

  if (tries > 1) {
    std::cout << Uintah::Parallel::getMPIRank()
              << " - WARNING:  Dir::createSubdirPlus() created the directory ("
              << sub << ") after " << tries << " attempts.\n";
  }
  return myDir;
}

Dir
Dir::getSubdir(const std::string& sub)
{
  // This should probably do more
  return Dir(d_name + "/" + sub);
}

void
Dir::copy_into(Dir& destDir)
{
  const auto copyOptions = fs::copy_options::recursive;
  auto last_dir          = fs::absolute(d_name).filename();
  auto dest_path         = fs::path(destDir.d_name);
  dest_path /= last_dir;
  try {
    fs::copy(fs::absolute(d_name), fs::absolute(dest_path), copyOptions);
  } catch (fs::filesystem_error& err) {
    throw InternalError(std::string("Dir::copy failed to copy: ") + d_name +
                          " to " + destDir.d_name + " " + err.what(),
                        __FILE__,
                        __LINE__);
  }
  return;
}

void
Dir::copy(Dir& destDir)
{
  const auto copyOptions = fs::copy_options::recursive;
  auto dest_path         = fs::path(destDir.d_name);
  try {
    fs::copy(fs::absolute(d_name), fs::absolute(dest_path), copyOptions);
  } catch (fs::filesystem_error& err) {
    throw InternalError(std::string("Dir::copy failed to copy: ") + d_name +
                          " to " + destDir.d_name + " " + err.what(),
                        __FILE__,
                        __LINE__);
  }
  return;
}

void
Dir::move(Dir& destDir)
{
  try {
    fs::rename(fs::absolute(d_name), fs::absolute(destDir.d_name));
  } catch (fs::filesystem_error& err) {
    auto last_dir  = fs::absolute(d_name).filename();
    auto dest_path = fs::path(destDir.d_name);
    dest_path /= last_dir;
    try {
      fs::rename(fs::absolute(d_name), fs::absolute(dest_path));
    } catch (fs::filesystem_error& err) {
      std::ostringstream msg;
      msg << "Dir::move failed to move: " << d_name << " to " << dest_path
          << " " << err.what();
      throw InternalError(msg.str(), __FILE__, __LINE__);
    }
  }
  return;
}

void
Dir::copy(const std::string& filename, Dir& destDir)
{
  fs::path filepath      = d_name;
  fs::path dest_filepath = destDir.d_name;

  filepath /= filename;
  try {
    dest_filepath /= filename;
    fs::copy_file(fs::absolute(filepath),
                  fs::absolute(dest_filepath),
                  fs::copy_options::overwrite_existing);
  } catch (fs::filesystem_error& err) {
    std::ostringstream msg;
    msg << "Dir::copy failed to copy file: " << fs::absolute(filepath) << " to "
        << fs::absolute(dest_filepath) << " with err " << err.what();
    throw InternalError(msg.str(), __FILE__, __LINE__);
  }
  return;
}

void
Dir::move(const std::string& filename, Dir& destDir)
{
  fs::path filepath      = d_name;
  fs::path dest_filepath = destDir.d_name;

  filepath /= filename;
  try {
    dest_filepath /= filename;
    fs::rename(fs::absolute(filepath), fs::absolute(dest_filepath));
  } catch (fs::filesystem_error& err) {
    std::ostringstream msg;
    msg << "Dir::move failed to move: " << fs::absolute(filepath) << " to "
        << fs::absolute(dest_filepath);
    throw InternalError(msg.str(), __FILE__, __LINE__);
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
Dir::exists() const
{
  if (fs::exists(d_name)) {
    return true;
  }
  return false;
}

} // namespace Uintah