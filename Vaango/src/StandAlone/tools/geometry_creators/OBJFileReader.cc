/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

//#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

// Stripped down OBJ file reader
//  Assumes ASCII file
//  Assumes first three vertex properties are x,y,z
//  Does not assume triangular faces - reads in number of vertices for each line

// Function declarations
void usage(const std::string& badarg, const std::string& progName);
std::string getFileWithoutExtension(const std::string& fileName);

int main(int argc, char** argv)
{
  // Parse arguments
  std::string input_file;
  std::string output_file_pts;
  std::string output_file_tri;

  // Read arguments and create output files
  if (argc < 2) {
    usage(" input file name not specified.", argv[0]);
  } else {
    std::cerr << "Using " << argv[1] << " as output file name" << std::endl;
    input_file = argv[1];
    std::string output_file = getFileWithoutExtension(input_file);
    output_file_pts = output_file + ".pts";
    output_file_tri = output_file + ".tri";
  }

  std::ifstream inputFile(argv[1]);
  if (!inputFile){
    std::cerr << "Error: Could not open input file " << argv[1] << std::endl;
    exit(1);
  }

  // Read the file
  std::string line;
  std::vector<double> xvert;
  std::vector<double> yvert;
  std::vector<double> zvert;
  std::vector<std::vector<int> > elements;

  while (std::getline(inputFile, line)) {

    // erase white spaces from the beginning of line
    line.erase(line.begin(), std::find_if(line.begin(), line.end(),
            [](unsigned char c){ return !std::isspace(c); }));

    // Ignore empty lines
    if (line.empty()) continue;

    // Skip comment lines
    if (line[0] == '#') continue;

    // Ignore lines containing texture information
    if (line[0] == 'v' && line[1] == 't') continue;

    // Read the vertices
    if (line[0] == 'v') {

      std::istringstream data_stream(line);
      char data_id;
      double xcoord, ycoord, zcoord = 0.0;
      if (!(data_stream >> data_id >> xcoord >> ycoord >> zcoord)) {
        std::cerr << "Could not read node input data stream" << std::endl;
        exit(1);
      }
      xvert.push_back(xcoord);
      yvert.push_back(ycoord);
      zvert.push_back(zcoord);
    }

    // Read the face connectivity
    if (line[0] == 'f') {
      std::istringstream data_stream4(line);
      std::istringstream data_stream3(line);
      char data_id;
      char dummy;
      int node1, node2, node3, node4;
      int tex1, tex2, tex3, tex4;
      std::vector<int> elem;
      if ((data_stream4 >> data_id >> node1 >> dummy >> tex1 
                                   >> node2 >> dummy >> tex2 
                                   >> node3 >> dummy >> tex3 
                                   >> node4 >> dummy >> tex4)) {
        std::cout << "[" << node1 << "," << node2 << "," << node3 << "," << node4 << "]" << std::endl;
        elem.push_back(node1);
        elem.push_back(node2);
        elem.push_back(node3);
        elem.push_back(node4);
        elements.push_back(elem);
        continue;
      } 
      if ((data_stream3 >> data_id >> node1 >> dummy >> tex1 
                                         >> node2 >> dummy >> tex2 
                                         >> node3 >> dummy >> tex3 )) {
        std::cout << "[" << node1 << "," << node2 << "," << node3 << "]" << std::endl;
        elem.push_back(node1);
        elem.push_back(node2);
        elem.push_back(node3);
        elements.push_back(elem);
        continue;
      } else {
        std::cout << data_stream3.str() << std::endl;
        std::cerr << "Could not read element input data stream" << std::endl;
        exit(1);
      }
    }
  }

  // Error if number of faces or vertices equals zero
  unsigned int n_vertex = xvert.size();
  unsigned int n_face = elements.size();
  if( n_vertex*n_face == 0 ) {
    std::cerr << "Error: Input file contains no vertices or faces" << std::endl;
    exit(1);
  }

  // Output vertices to .pts file
  std::cout << "Writing " << n_vertex << " vertices" << std::endl;
  std::ofstream pts_stream;
  pts_stream.open(output_file_pts.c_str());
  for(unsigned int ii=0; ii < n_vertex; ii++ ) {
    pts_stream << xvert[ii] << " " << yvert[ii] << " " << zvert[ii] << std::endl;
  }
  pts_stream.close();

  // Output faces to .tri file
  std::cout << "Writing " << n_face << " faces" << std::endl;
  std::ofstream tri_stream;
  tri_stream.open(output_file_tri.c_str());
  for (unsigned int ii=0; ii< n_face; ii++ ) {
    if (elements[ii].size() == 3) {
      for (auto iter = elements[ii].begin(); iter != elements[ii].end(); iter++) {
        tri_stream << *iter << " ";
      }
      tri_stream << std::endl;
    } else if (elements[ii].size() == 4) {
      std::vector<int> elem = elements[ii];
      tri_stream << elem[0] << " " << elem[1] << " " << elem[2] << std::endl;
      tri_stream << elem[0] << " " << elem[2] << " " << elem[3] << std::endl;
    } else {
      std::cerr << "Error: Input file contains faces with < 3 or > 4 faces" << std::endl;
      exit(1);
    }
  }
  tri_stream.close();

  return (0);
}

// Usage
void usage(const std::string& badarg, const std::string& progname)
{
  if(badarg != "") std::cerr << "Error parsing argument: " << badarg << std::endl;

  std::cerr << "Usage: " << progname << " <input file>\n\n";
  std::cerr << "  e.g., OBJFileReader sphere.obj " << std::endl;
  std::cerr << "        will create two files: sphere.pts and sphere.tri" << std::endl;
  exit(1);
}

std::string getFileWithoutExtension(const std::string& fileName)
{
  if (fileName.find_last_of(".") != std::string::npos)
    return fileName.substr(0, fileName.find_last_of("."));
  return "";
}
