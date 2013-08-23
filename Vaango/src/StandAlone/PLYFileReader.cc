#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

// Stripped down PLY reader
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
	if(!inputFile){
		std::cerr << "Error: Could not open input file" << std::endl;
		exit(1);
	}

  // Read header - get number of vertices and faces
  std::string in_string;
  std::string end_head = "end_header";
  int n_vertex, n_face;
  do {
    inputFile >> in_string;
    if( !in_string.compare("vertex") ) inputFile >> n_vertex;
    if( !in_string.compare("face") ) inputFile >> n_face;
	} while (!inputFile.eof() && in_string.compare(end_head) );

  // Error if number of faces or vertices equals zero
  if( n_vertex*n_face == 0 ) {
    std::cerr << "Error: Input file contains no vertices or faces" << std::endl;
    exit(1);
  }

  // Read in each vertex - output to .pts file
  std::cout << "Loading " << n_vertex << " vertices" << std::endl;
  double x, y, z;
  std::string str_tmp;
  std::ofstream pts_stream;
  pts_stream.open(output_file_pts.c_str());
  for( int ii=0; ii<n_vertex; ii++ ) {
    inputFile >> x >> y >> z;
    std::getline( inputFile, str_tmp );
    pts_stream << x << " " << y << " " << z << std::endl;
  }
  pts_stream.close();

  // Read in each face - output to .tri file
  std::cout << "Loading " << n_face << " faces" << std::endl;
  int n_vtx, vtx;
  std::ofstream tri_stream;
  tri_stream.open(output_file_tri.c_str());
  for( int ii=0; ii<n_face; ii++ ) {
    inputFile >> n_vtx;
    for( int ivtx=0; ivtx<n_vtx; ivtx++ ) {
      inputFile >> vtx;
      tri_stream << vtx << " ";
    }
    tri_stream << std::endl;
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
