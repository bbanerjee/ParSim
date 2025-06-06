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

// Test triangulation of a surface
#include <Core/Exception.h>

#include <pcl/surface/organized_fast_mesh.h>
#include <pcl/point_types.h>
#include <pcl/io/vtk_io.h>

#include <chrono>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using namespace Matiti;

void test_organized();

int main()
{
  auto t1 = std::chrono::high_resolution_clock::now();
  test_organized();
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Time: " 
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
            << std::endl;
  return 0;
}

void test_organized()
{
  typedef pcl::PointCloud<pcl::PointXYZ> Cloud;
  typedef Cloud::Ptr                     CloudP;
  
  // create empty point cloud
  CloudP cloud(new Cloud());

  // Read in the points
  std::string fileName("bunny_sorted.dat");

  // Try to open file
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::string out = "Could not open point input file " + fileName + " for reading \n";
    throw Exception(out, __FILE__, __LINE__);
  }

  // Read the file and fille the cloud
  std::string line;
  while (std::getline(file, line)) {

    // erase white spaces from the beginning of line
    line.erase(line.begin(), std::find_if(line.begin(), line.end(),
         std::not1(std::ptr_fun<int, int>(std::isspace))));

    // Read the data
    std::istringstream data_stream(line);
    double xcoord, ycoord;
    if (!(data_stream >> xcoord >> ycoord)) {
      throw Exception("Could not read point input data stream", __FILE__, __LINE__);
    }
    cloud->push_back(pcl::PointXYZ(xcoord, ycoord, 0.0));
  }

  // Set up meshing objects
  pcl::PolygonMesh triangles;
  pcl::OrganizedFastMesh<pcl::PointXYZ> mesh;

  // Initialize
  mesh.setInputCloud(cloud);
  mesh.setMaxEdgeLength(10);
  mesh.setTriangulationType(pcl::OrganizedFastMesh<pcl::PointXYZ>::TRIANGLE_ADAPTIVE_CUT);

  // Reconstruct
  mesh.reconstruct(triangles);

  // Save triangles
  pcl::io::saveVTKFile("./test.vtk", triangles);

}

