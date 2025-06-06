/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 1997-2012 The University of Utah
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

/* FileUtils.cc */
#include <Core/Util/FileUtils.h>

#include <Core/OS/Dir.h>
#include <Core/Util/Assert.h>
#include <Core/Util/Environment.h>
#include <Core/Util/StringUtil.h>

#include <iostream>
#include <sstream>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>

namespace Uintah {

////////////////////////////////////////////////////////////
//
// InsertStringInFile
//
//       If "match" is found in "filename", then add "add_text" text
//       to the file _after_ the location of the matched text.
//
//   Normally, I would just use sed via system() to edit a file,
//   but for some reason system() calls never work from Dataflow
//   processes in linux.  Oh well, sed isn't natively available
//   under windows, so I'd have to do something like this anyhow
//   - Chris Moulding

void
InsertStringInFile(char* filename, const char* match, const char* add_text)
{
  char* newfilename = new char[strlen(filename) + 2];
  char c;
  sprintf(newfilename, "%s~", filename);
  FILE* ifile;
  FILE* ofile;

  /* create a copy of the original file */
  ifile = fopen(filename, "r");

  if (ifile == nullptr) {
    printf("ERROR: In Core/Util/FileUtils.cc: InsertStringInFile:\n");
    printf("        File '%s' does not exist!\n", filename);
    printf("       There is something seriously wrong with your Uintah "
           "installation.\n");
    printf("       Please contact scirun@sci.utah.edu.\n");
    exit(1);
  }

  ofile = fopen(newfilename, "w");

  c = (char)fgetc(ifile);
  while (c != (char)EOF) {
    fprintf(ofile, "%c", c);
    c = (char)fgetc(ifile);
  }
  fclose(ifile);
  fclose(ofile);

  // Search the copy for an instance of "match"...
  int index1          = 0;
  unsigned int index2 = 0;
  int foundat         = -1;
  ifile               = fopen(newfilename, "r");
  c                   = (char)fgetc(ifile);
  while (c != (char)EOF) {
    if (c == match[index2]) {
      foundat = index1 + strlen(match);
      while (index2 < strlen(match) && c != (char)EOF && c == match[index2]) {
        c = (char)fgetc(ifile);
        index1++;
        index2++;
      }
      if (foundat >= 0 && index2 != strlen(match)) {
        foundat = -1;
        index2  = 0;
      } else {
        break;
      }
    }
    c = (char)fgetc(ifile);
    index1++;
  }
  fclose(ifile);

  // If an instance of match was found, insert the indicated string...
  if (foundat >= 0) {
    index1 = 0;
    ifile  = fopen(newfilename, "r");
    ofile  = fopen(filename, "w");
    c      = (char)fgetc(ifile);
    while (c != (char)EOF) {
      if (index1 == foundat) {
        fprintf(ofile, "%s", add_text);
      }
      fprintf(ofile, "%c", c);
      c = (char)fgetc(ifile);
      index1++;
    }
    fclose(ifile);
    fclose(ofile);
  }
}

std::map<int, char*>*
GetFilenamesEndingWith(const char* d, std::string ext)
{
  std::map<int, char*>* newmap = 0;
  dirent* file                 = 0;
  DIR* dir                     = opendir(d);
  char* newstring              = 0;

  if (!dir) {
    return 0;
  }

  newmap = new std::map<int, char*>;

  file = readdir(dir);
  while (file) {
    if ((strlen(file->d_name) >= strlen(ext.c_str())) &&
        (strcmp(&(file->d_name[strlen(file->d_name) - strlen(ext.c_str())]),
                ext.c_str()) == 0)) {
      newstring = new char[strlen(file->d_name) + 1];
      sprintf(newstring, "%s", file->d_name);
      newmap->insert(std::pair<int, char*>(newmap->size(), newstring));
    }
    file = readdir(dir);
  }

  closedir(dir);
  return newmap;
}

std::string
substituteTilde(const std::string& dirstr)
{
  std::string realdirstr = dirstr;

  std::string::size_type pos = realdirstr.find("~");
  if (pos != std::string::npos && (pos == 0 || dirstr[pos - 1] == '/')) {
    std::string HOME = sci_getenv("HOME");
    realdirstr = HOME + "/" + dirstr.substr(pos + 1, dirstr.size() - pos - 1);
  }
  return realdirstr;
}

std::vector<std::string>
GetFilenamesStartingWith(const std::string& dirstr, const std::string& prefix)
{
  std::vector<std::string> files(0);
  DIR* dir = opendir(substituteTilde(dirstr).c_str());
  if (!dir) {
    return files;
  }

  dirent* file = readdir(dir);
  while (file) {
    ASSERT(file->d_name);
    std::string dname          = file->d_name;
    std::string::size_type pos = dname.find(prefix);
    if (pos == 0) {
      files.push_back(dname);
    }
    file = readdir(dir);
  }

  closedir(dir);
  return files;
}

std::pair<std::string, std::string>
split_filename(const std::string& fname)
{
  std::string filename = fname;

  if (fname[fname.size() - 1] == '/' || fname[fname.size() - 1] == '\\') {
    filename = fname.substr(0, fname.size() - 1);
  }

  if (validDir(fname)) {
    return make_pair(filename, std::string(""));
  }

  std::string::size_type pos = filename.find_last_of("/");
  if (pos == std::string::npos) {
    pos = filename.find_last_of("\\");
  }
  std::pair<std::string, std::string> dirfile =
    make_pair(filename.substr(0, pos + 1),
              filename.substr(pos + 1, filename.length() - pos - 1));

  return dirfile;
}

bool
getInfo(const std::string& filename)
{
  struct stat buf;
  std::string updatedFilename = substituteTilde(filename);

  int result = stat(updatedFilename.c_str(), &buf);
  if (result == 0) {
    printf("Successful stat of file '%s'.\n", filename.c_str());
    mode_t& m = buf.st_mode;

    if (m & S_IRUSR && S_ISREG(m) && !S_ISDIR(m)) {
      printf("   File appears to be a regular file.\n");
    }
    if (m & S_IRUSR && !S_ISREG(m) && S_ISDIR(m)) {
      printf("   File appears to be a regular directory.\n");
    }
    if (m & S_ISLNK(m)) {
      printf("   File appears to be a symbolic link.\n");
    }
    printf("   File size: %d\n", (int)buf.st_size);

    struct tm* tmptr;
    tmptr = localtime(&buf.st_mtime);
    char time_str[256];
    strftime(time_str, 250, "%D %T", tmptr);

    printf("   Last file modification: %s\n", time_str);
    return true;
  } else {
    printf("Error getting information on file '%s'.  Errno: %d.\n",
           filename.c_str(),
           errno);
    perror("   Reason");
    return false;
  }
}

bool
validFile(const std::string& filename)
{
  struct stat buf;
  std::string updatedFilename = substituteTilde(filename);

  int result = stat(updatedFilename.c_str(), &buf);
  if (result == 0) {
    mode_t& m  = buf.st_mode;
    bool valid = m & S_IRUSR && S_ISREG(m) && !S_ISDIR(m);
    return valid;
  }
  return false;
}

bool
validDir(const std::string& dirname)
{
  struct stat buf;
  std::string updatedDirname = substituteTilde(dirname);
  if (stat(updatedDirname.c_str(), &buf) == 0) {
    mode_t& m = buf.st_mode;
    return (m & S_IRUSR && !S_ISREG(m) && S_ISDIR(m));
  }
  return false;
}

bool
isSymLink(const std::string& filename)
{
  struct stat buf;
  if (lstat(filename.c_str(), &buf) == 0) {
    mode_t& m = buf.st_mode;
    return (m & S_ISLNK(m));
  }
  return false;
}

// Creates a temp file (in directoryPath), writes to it, and then deletes it...
bool
testFilesystem(const std::string& directoryPath,
               std::stringstream& error_stream,
               int procNumber /* = -1 ... ie, non-MPI test...*/)
{
  FILE* fp;

  std::string fileName = directoryPath + "/scirun_filesystem_check_temp_file";
  if (procNumber != -1) {
    std::stringstream fn;
    fn << fileName;
    fn << "_" << procNumber;
    fileName = fn.str();
  }

  // Create a temporary file
  fp = fopen(fileName.c_str(), "w");
  if (fp == nullptr) {
    error_stream << "ERROR: testFilesystem() failed to create a temp file ("
                 << fileName << ") in " << directoryPath << "\n";
    error_stream << "       errno is " << errno << "\n";
    return false;
  }

  // Write to the file
  const char* myStr = "hello world";
  for (int cnt = 0; cnt < 1001; cnt++) {
    int numWritten = fwrite(myStr, 1, 11, fp);
    if (numWritten != 11) {
      error_stream
        << "ERROR: testFilesystem() failed to write data to temp file ("
        << fileName << ") in " << directoryPath << "\n";
      error_stream << "       iteration: " << cnt << ", errno is " << errno
                   << "\n";
      return false;
    }
  }

  // Close the file
  int result = fclose(fp);
  if (result != 0) {
    error_stream << "WARNING: fclose() failed while testing filesystem.\n";
    error_stream << "         errno is " << errno << "\n";
    return false;
  }

  // Check the files size
  struct stat buf;
  if (stat(fileName.c_str(), &buf) == 0) {
    if (buf.st_size != 11011) {
      error_stream << "ERROR: Test file size is: " << (int)buf.st_size
                   << ", but should be 11011 bytes!\n";
    }
  } else {
    error_stream << "ERROR: stat() failed while testing filesystem.\n";
    error_stream << "         errno is " << errno << "\n";
    return false;
  }

  // Delete the file
  int rc = remove(fileName.c_str());
  if (rc != 0) {
    error_stream << "ERROR: remove() failed while testing filesystem.\n";
    error_stream << "         errno is " << errno << "\n";
    return false;
  }

  return true;
}

// findFileInPath
// Searches the colon-seperated 'path' std::string variable for a file named
// 'file'.  From left to right in the path string each directory is
// tested to see if the file named 'file' is in it.
//
// If the file is found, it returns the DIRECTORY that the file is located in
// Otherwise if the file is not found in the path, returns an empty string
//
// If the file is found in multiple directories in the 'path', only
// the first matching directory is returned
std::string
findFileInPath(const std::string& file, const std::string& path)
{
  std::string::size_type beg = 0;
  std::string::size_type end = 0;

  while (beg < path.length()) {
    end = path.find(":", beg);
    if (end == std::string::npos) {
      end = path.size();
    }

    std::string dir(path, beg, end - beg);
    ASSERT(!dir.empty());
    // Append a slash if there isn't one
    if (validDir(dir)) {
      if (dir[dir.length() - 1] != '/') {
        dir = dir + "/";
      }
      std::string filename = dir + file;
      if (validFile(filename)) {
        // see comments at function start
        return dir;
      }
    }

    beg = end + 1;
  }
  return "";
}

std::string
autocomplete(const std::string& instr)
{
  std::string str                             = instr;
  std::pair<std::string, std::string> dirfile = split_filename(str);
  std::string dir                             = dirfile.first;
  if (!validDir(dir)) {
    return str;
  }
  if (dirfile.first.empty() ||
      dirfile.first[dirfile.first.length() - 1] != '/') {
    dirfile.first = dirfile.first + "/";
  }
  std::vector<std::string> files = GetFilenamesStartingWith(dir, dirfile.second);

  if (files.empty()) {
    return str;
  }
  if (files.size() == 3 && files[0] == "." && files[1] == "..") {
    str = dirfile.first + files[2];
    if (validDir(str)) {
      str = str + "/";
    }
  } else if (files.size() == 1) {
    str = dirfile.first + files[0];
    if (validDir(str)) {
      str = str + "/";
    }
  } else {
    unsigned int j0 = dirfile.second.size();
    unsigned int j  = j0;
    do {
      for (unsigned int i = 1; i < files.size(); ++i) {
        if ((j == files[i].size()) || (files[i][j] != files[i - 1][j])) {
          str = str + files[i].substr(j0, j - j0);
          return str;
        }
      }
      ++j;
    } while (1);
  }
  return str;
}

std::string
canonicalize(const std::string& fname)
{
  std::string filename = substituteTilde(fname);

  // use unix paths internally to keep things simpler
  convertToUnixPath(filename);
  replace_substring(filename, "//", "/");

  std::vector<char> separators;
  separators.push_back('\'');

  std::vector<std::string> entries, newentries = split_string(filename, separators);
  while (newentries.size() != entries.size()) {
    entries = newentries;
    newentries.clear();
    for (unsigned int i = 0; i < entries.size(); ++i) {
      if ((entries[i] != ".") &&  // Skip . entries
          (entries[i] != "..") && // Skil .. entries
          (i + 1 == entries.size() ||
           entries[i + 1] != "..")) // Skip entries parenting .. entries
      {
        newentries.push_back(entries[i]);
      }
    }
  }

  filename = "";
  for (unsigned int i = 0; i < entries.size(); ++i) {
    filename = filename + "/" + entries[i];
  }
  if (filename == "" || validDir(filename)) {
    filename = filename + "/";
  }

  return filename;
}

void
convertToWindowsPath(std::string& unixPath)
{
  for (std::string::size_type cnt = 0; cnt < unixPath.length(); cnt++) {
    if (unixPath[cnt] == '/') {
      unixPath[cnt] = '\\';
    }
  }
}

void
convertToUnixPath(std::string& unixPath)
{
  for (std::string::size_type cnt = 0; cnt < unixPath.length(); cnt++) {
    if (unixPath[cnt] == '\\') {
      unixPath[cnt] = '/';
    }
  }
}

int
copyFile(const std::string& src, const std::string& dest)
{
  std::string cpCmd = "cp -f ";
  std::string cmd   = cpCmd + src + " " + dest;
  int code          = std::system(cmd.c_str());
  if (code) {
    std::cerr << "Error executing: " << cmd << "\n";
  }
  return code;
}

int
moveFile(const std::string& src, const std::string& dest)
{
  std::string mvCmd = "mv -f ";
  std::string cmd   = mvCmd + src + " " + dest;
  int code          = std::system(cmd.c_str());
  if (code) {
    std::cerr << "Error executing: " << cmd << "\n";
  }
  return code;
}

int
deleteFile(const std::string& filename)
{
  std::string rmCmd = "rm -f ";
  std::string cmd   = rmCmd + filename;
  int code          = std::system(cmd.c_str());
  if (code) {
    std::cerr << "Error executing: " << cmd << "\n";
  }
  return code;
}

int
copyDir(const std::string& src, const std::string& dest)
{
  std::string cpCmd = "cp -fr ";
  std::string cmd   = cpCmd + src + " " + dest;
  int code          = std::system(cmd.c_str());
  if (code) {
    std::cerr << "Error executing: " << cmd << "\n";
  }
  return code;
}

int
deleteDir(std::string filename)
{
  std::string rmCmd = "rm -rf ";
  std::string cmd   = rmCmd + filename;
  int code          = std::system(cmd.c_str());
  if (code) {
    std::cerr << "Error executing: " << cmd << "\n";
  }
  return code;
}

std::string
changeExtension(const std::string& filename, const std::string& extension)
{
  std::string fname = filename;

  std::pair<std::string, std::string> dirfile = split_filename(filename);
  if (dirfile.second.size()) {

    std::vector<char> separators;
    separators.push_back('.');

    std::vector<std::string> fileext = split_string(dirfile.second, separators);
    fname                       = dirfile.first;
    for (size_t i = 0; i < fileext.size() - 1; ++i) {
      fname = fname + fileext[i] + ".";
    }
    fname = fname + extension;
  }
  return fname;
}

} // End namespace Uintah
