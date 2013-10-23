#include <GeometryPiece/BoxGeometryPiece.h>
#include <Exception.h>
#include <MPMDatawarehouse.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Geometry/Vector.h>

#include <iostream>

using namespace BrMPM; 

BoxGeometryPiece::BoxGeometryPiece(Uintah::ProblemSpecP& ps,
                                   MPMDatawarehouseP& dw)
{
  d_name = "box";
  Uintah::Vector lower, upper;
  ps->require("lower", lower);
  ps->require("upper", upper);
  d_box = Box3D(Point3D(lower[0],lower[1],lower[2]), Point3D(upper[0],upper[1],upper[2]));
  if (d_box.isDegenerate()) {
    std::ostringstream out;
    out << "**ERROR** The box geometry piece is degenerate. Lower = " << d_box.lower()
        << " Upper = " << d_box.upper() << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  // Create nodes
  createParticles(dw);
}

BoxGeometryPiece::~BoxGeometryPiece()
{
}

Box3D 
BoxGeometryPiece::boundingBox() const
{
  return d_box;
}


bool 
BoxGeometryPiece::inside (const Point3D& pt) const
{
  return d_box.contains(pt); 
}


std::string 
BoxGeometryPiece::name() const
{
  return d_name;
}

void
BoxGeometryPiece::createParticles(MPMDatawarehouseP& dw)
{
  // TODO: Hooman please complete this

  /*
  // Create nodes
  int nx = d_num_elements[0];
  int ny = d_num_elements[1];
  int nz = d_num_elements[2];
  double xmin = (d_box.lower())[0];
  double ymin = (d_box.lower())[1];
  double zmin = (d_box.lower())[2];
  Vector3D span = d_box.upper() - d_box.lower();
  double dx = span.x()/(double) nx;
  double dy = span.y()/(double) ny;
  double dz = span.z()/(double) nz;

  std::vector<double> xcoords, ycoords, zcoords;
  nx++;
  for (int ii=0; ii < nx; ++ii) {
    xcoords.emplace_back(xmin + dx*(double) ii);
  }
  ny++;
  for (int jj=0; jj < ny; ++jj) {
    ycoords.emplace_back(ymin + dy*(double) jj);
  }
  nz++;
  for (int kk=0; kk < nz; ++kk) {
    zcoords.emplace_back(zmin + dz*(double) kk);
  }
  int node_id = 0;
  int node_z = 0;
  for (auto ziter = zcoords.begin(); ziter != zcoords.end(); ++ziter) {
    ++node_z;
    bool on_z_surf = (node_z == 1) || (node_z == nz);
    int node_y = 0;
    for (auto yiter = ycoords.begin(); yiter != ycoords.end(); ++yiter) {
      ++node_y;
      bool on_y_surf = (node_y == 1) || (node_y == ny);
      int node_x = 0;
      for (auto xiter = xcoords.begin(); xiter != xcoords.end(); ++xiter) {
        ++node_x;
        bool on_x_surf = (node_x == 1) || (node_x == nx);
        ++node_id;
        NodeP node(new Node(node_id, *xiter, *yiter, *ziter, on_x_surf||on_y_surf||on_z_surf));
        nodes.emplace_back(node);
        //std::cout << "node = " << node_id << " " <<  node.position() ;

       // Add to the node ID -> node ptr map
       d_id_ptr_map.insert(std::pair<int, NodeP>(node_id, node));
      }
    }
  }
  */
}

