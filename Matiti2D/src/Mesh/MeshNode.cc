#include <Vaango/Core/Geometry/Point.h>
#include <Mesh/MeshNode.h>
#include <Mesh/Point.h>
#include <iostream>

using namespace std;
using namespace Matiti;

MeshNode::MeshNode(const long int& id, const double& x, const double& y)
{
	d_id = id;
	d_x = x;
	d_y = y;
	d_z = 0.0;
}

MeshNode::MeshNode(const long int& id, const double& x, const double& y, const double& z)
{
	d_id = id;
	d_x = x;
	d_y = y;
	d_z = z;
}

MeshNode::MeshNode(const long int& id, const SCIRun::Point& pt)
{
	d_id = id;
	d_x = pt.x_;
	d_y = pt.y_;
	d_z = pt.z_;
}

MeshNode:: ~MeshNode()
{
}

long int MeshNode::id() const
{
	return d_id;
}

SCIRun::Point& MeshNode::position() const
{
  return SCIRun::Point(d_x, d_y, d_z);
}
