/*
 * The MIT License
 *
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

#include "happly.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <regex>

// Stripped down PLY reader
//  Assumes ASCII file
//  Assumes first three vertex properties are x,y,z
//  Does not assume triangular faces - reads in number of vertices for each line

// Function declarations
void usage(const std::string& badarg, const std::string& progName);
std::string getFileWithoutExtension(const std::string& fileName);
void readBinaryPLYFile(const std::string& fileName, const std::string& tri, const std::string& pts);

int
main(int argc, char** argv)
{
  // Parse arguments
  std::string input_file;
  std::string output_file_pts;
  std::string output_file_tri;

  // Read arguments and create output files
  if (argc < 2) {
    usage(" input file name not specified.", argv[0]);
  } else {
    std::cerr << "Using " << argv[1] << " as output file name" << "\n";
    input_file = argv[1];
    std::string output_file = getFileWithoutExtension(input_file);
    output_file_pts = output_file + ".pts";
    output_file_tri = output_file + ".tri";
  }

  std::ifstream inputFile(argv[1]);
  if (!inputFile) {
    std::cerr << "Error: Could not open input file" << "\n";
    exit(1);
  }

  // Check if this is a PLY file
  // Check if the format is binary or ascii
  std::string in_string;
  const std::regex one_or_more_space("\\s+");
  bool plyFile = false; 
  bool asciiFile = false; 
  int count = 0;
  while (std::getline(inputFile, in_string)) {
    if (in_string.empty()) continue;
    std::sregex_token_iterator it(in_string.begin(), in_string.end(), one_or_more_space, -1); 
    if (*it == "ply") {
      plyFile = true;
    }
    if (*it == "format") {
      ++it;
      if (*it == "ascii") {
        asciiFile = true;
      }
    } 
    ++count;
    if (count > 2) break;
  }
  if (!plyFile) {
    inputFile.close();
    std::cerr << "Error: Input file is not a PLY file." << "\n";
    exit(1);
  }
  if (!asciiFile) {
    inputFile.close();
    std::cerr << "Error: Input file is not an ASCII PLY file." << "\n";
    readBinaryPLYFile(input_file, output_file_pts, output_file_tri);
    exit(1);
  }
  
  // Read header - get number of vertices and faces
  std::string end_head = "end_header";
  int n_vertex, n_face;
  do {
    inputFile >> in_string;
    if (!in_string.compare("vertex"))
      inputFile >> n_vertex;
    if (!in_string.compare("face"))
      inputFile >> n_face;
  } while (!inputFile.eof() && in_string.compare(end_head));

  // Error if number of faces or vertices equals zero
  if (n_vertex * n_face == 0) {
    std::cerr << "Error: Input file contains no vertices or faces" << "\n";
    exit(1);
  }

  // Read in each vertex - output to .pts file
  std::cout << "Loading " << n_vertex << " vertices" << "\n";
  double x, y, z;
  std::string str_tmp;
  std::ofstream pts_stream;
  pts_stream.open(output_file_pts.c_str());
  for (int ii = 0; ii < n_vertex; ii++) {
    inputFile >> x >> y >> z;
    std::getline(inputFile, str_tmp);
    pts_stream << x << " " << y << " " << z << "\n";
  }
  pts_stream.close();

  // Read in each face - output to .tri file
  std::cout << "Loading " << n_face << " faces" << "\n";
  int n_vtx, vtx;
  std::ofstream tri_stream;
  tri_stream.open(output_file_tri.c_str());
  for (int ii = 0; ii < n_face; ii++) {
    inputFile >> n_vtx;
    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      inputFile >> vtx;
      tri_stream << vtx << " ";
    }
    tri_stream << "\n";
  }
  tri_stream.close();

  inputFile.close();

  return (0);
}

void 
readBinaryPLYFile(const std::string& inFile,
                  const std::string& outFilePts,
                  const std::string& outFileTri)
{
  // Read the file
  happly::PLYData plyIn(inFile);
  std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
  std::vector<std::vector<size_t>> faces = plyIn.getFaceIndices<size_t>();

  // Output to .pts file
  std::cout << "Writing " << outFilePts << "\n";
  std::ofstream pts_stream;
  pts_stream.open(outFilePts.c_str());
  for (auto pos : vPos) {
    pts_stream << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
  }
  pts_stream.close();

  // Output to .tri file
  std::cout << "Writing " << outFileTri << "\n";
  std::ofstream tri_stream;
  tri_stream.open(outFileTri.c_str());
  for (auto face : faces) {
    for (auto index : face) {
      tri_stream << index << " " ;
    }
    tri_stream << "\n";
  }
  tri_stream.close();
}

// Usage
void
usage(const std::string& badarg, const std::string& progname)
{
  if (badarg != "")
    std::cerr << "Error parsing argument: " << badarg << "\n";

  std::cerr << "Usage: " << progname << " <input file>\n\n";
  std::cerr << "  e.g., PLYFileReader bunny.ply " << "\n";
  std::cerr << "        will create two files: bunny.pts and bunny.tri"
            << "\n";
  exit(1);
}

std::string
getFileWithoutExtension(const std::string& fileName)
{
  if (fileName.find_last_of(".") != std::string::npos)
    return fileName.substr(0, fileName.find_last_of("."));
  return "";
}
