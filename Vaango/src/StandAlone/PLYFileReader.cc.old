#include <iostream>
#include <fstream>
#include <string>
#include <pcl/io/ply_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/point_types.h>

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

  // Read point cloud from PLY file
  sensor_msgs::PointCloud2 cloud2;
  //if (reader.read(input_file, cloud2) < 0) {
  if (pcl::io::loadPLYFile(input_file, cloud2) < 0) {
    std::cout << "Couldn't read " << input_file << std::endl;
  }

  // Convert from PointCloud2 to pcl::PointColud<T>
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
  pcl::fromROSMsg(cloud2, *cloud);

  // Print the point data
  ofstream pts_stream;
  pts_stream.open(output_file_pts.c_str());
  std::cout << "Loaded "
            << cloud->width * cloud->height
            << " data points from test_pcd.pcd with the following fields: "
            << std::endl;
  for (size_t i = 0; i < cloud->points.size (); ++i) {
    //std::cout << "    " << cloud->points[i].x
    //          << " "    << cloud->points[i].y
    //          << " "    << cloud->points[i].z << std::endl;
    pts_stream << cloud->points[i].x << " " << cloud->points[i].y
               << " " << cloud->points[i].z << std::endl;
  }
  pts_stream.close();


  // Read polygon mesh from PLY file
  pcl::PolygonMesh mesh;
  if (pcl::io::loadPolygonFilePLY(input_file, mesh) < 0) {
    std::cout << "Couldn't read " << input_file << std::endl;
  }

  // Print the polygon data
  ofstream tri_stream;
  tri_stream.open(output_file_tri.c_str());
  std::cout << "Loaded "
            << mesh.cloud.width * mesh.cloud.height
            << " data points from test_pcd.pcd with the following fields: "
            << std::endl;
  for (size_t i = 0; i < mesh.polygons.size (); ++i) {
    //std::cout << "    " << mesh.polygons[i].vertices[0] 
    //          << " "    << mesh.polygons[i].vertices[1] 
    //          << " "    << mesh.polygons[i].vertices[2] << std::endl;
    tri_stream << mesh.polygons[i].vertices[0] << " "    
               << mesh.polygons[i].vertices[1] << " "    
               << mesh.polygons[i].vertices[2] << std::endl;
  }
  tri_stream.close();

  return (0);
  
}

// Usage
void usage(const std::string& badarg, const std::string& progname)
{
  if(badarg != "") std::cerr << "Error parsing argument: " << badarg << std::endl;

  std::cerr << "Usage: " << progname << " <input file>\n\n";
  std::cerr << "  e.g., PLYFileReader bunny.ply " << std::endl;
  std::cerr << "        will create two files: bunny.pts and bunny.tri" << std::endl;
  exit(1);
}

std::string getFileWithoutExtension(const std::string& fileName)
{
  if (fileName.find_last_of(".") != std::string::npos)
    return fileName.substr(0, fileName.find_last_of("."));
  return "";
}
