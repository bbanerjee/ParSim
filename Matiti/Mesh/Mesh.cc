#include <Core/Domain/Mesh.h>
#include <Core/Domain/Level.h>
#include <Core/Domain/Domain.h>
#include <Core/Domain/Variables/MeshElementIterator.h>
#include <Core/Domain/Variables/MeshNodeIterator.h>
#include <Core/Exceptions/InvalidDomain.h>
#include <Core/Math/Primes.h>
#include <Core/Domain/Box.h>
#include <Core/Domain/BoundaryConditions/BCData.h>
#include <Core/Domain/BoundaryConditions/BCDataArray.h>
#include <Core/Domain/BoundaryConditions/BoundCond.h>
#include <Core/Containers/StaticArray.h>
#include <TauProfilerForSCIRun.h>
#include <Core/Thread/AtomicCounter.h>
#include <Core/Thread/Mutex.h>
#include <Core/Math/MiscMath.h>

#include <iostream>
#include <sstream>
#include <cstdio>
#include <map>

using namespace std;
using namespace SCIRun;
using namespace Matiti;


static AtomicCounter ids("Mesh ID counter",0);
static Mutex ids_init("ID init");


Mesh::Mesh(int id):d_domain(0)
{
  
  if(d_id == -1){
    d_id = ids++;

  } else {
    if(d_id >= ids)
      ids.set(d_id+1);
  }
}

Mesh::Mesh(const Mesh* realMesh)
    : d_domain(realMesh->d_domain),
{
  d_id = realMesh->d_id; 
}

Mesh::~Mesh()
{
}

/**
* Returns the 8 meshNodes found around the point pos
*/
void Mesh::findMeshElementMeshNodes(const Point& pos, IntVector ni[8]) const
{
}

/**
 * Returns the 27 meshNodes found around the point pos
 */
void Mesh::findMeshElementMeshNodes27(const Point& pos, IntVector ni[27]) const
{
}


/**
 * Returns the position of the meshNode idx in domain coordinates.
 */
Point Mesh::meshNodePosition(const IntVector& idx) const {
}

/**
 * Returns the position of the meshElement idx in domain coordinates.
 */
Point Mesh::meshElementPosition(const IntVector& idx) const {
}

void Mesh::findMeshElementsFromMeshNode( const IntVector& meshNodeIndex,
                               IntVector meshElementIndex[8]) 
{
}

void Mesh::findMeshNodesFromMeshElement( const IntVector& meshElementIndex,
                               IntVector meshNodeIndex[8])
{
}

namespace Matiti {
  ostream&
  operator<<(ostream& out, const Mesh & r)
  {
    out.setf(ios::scientific,ios::floatfield);
    out.precision(4);
    out.setf(ios::scientific ,ios::floatfield);
    return out;
  }
}

//void
//Mesh::setBCType(Mesh::FaceType face, BCType newbc)
//{
//}


/**
 * This function will return an iterator that touches all meshElements
 *  that are are partially or fully within the region formed
 *  by intersecting the box b and this mesh (including extra meshElements).
 *  The region is inclusive on the + faces
 */
MeshElementIterator
Mesh::getMeshElementIterator(const Box& b) const
{
   Point l = getLevel()->positionToIndex(b.lower());
   Point u = getLevel()->positionToIndex(b.upper());
   IntVector low(RoundDown(l.x()), RoundDown(l.y()), RoundDown(l.z()));
   // high is the inclusive upper bound on the index.  In order for
   // the iterator to work properly we need in increment all the
   // indices by 1.
   IntVector high(RoundDown(u.x())+1, RoundDown(u.y())+1, RoundDown(u.z())+1);
   low = Max(low, getExtraMeshElementLowIndex());
   high = Min(high, getExtraMeshElementHighIndex());
   return MeshElementIterator(low, high);
}

/**
 *  This function will return an iterator that touches all meshElements
 *  whose centers are within the region formed
 *  by intersecting the box b and this mesh (including extra meshElements).
 *  The region is inclusive on the + faces
 */
MeshElementIterator
Mesh::getMeshElementCenterIterator(const Box& b) const
{
  Point l = getLevel()->positionToIndex(b.lower());
  Point u = getLevel()->positionToIndex(b.upper());
  // If we subtract 0.5 from the bounding box locations we can treat
  // the code just like we treat meshNodes.
  l -= Vector(0.5, 0.5, 0.5);
  u -= Vector(0.5, 0.5, 0.5);
  // This will return an empty iterator when the box is degerate.
  IntVector low(RoundUp(l.x()), RoundUp(l.y()), RoundUp(l.z()));
  IntVector high(RoundDown(u.x()) + 1, RoundDown(u.y()) + 1,
      RoundDown(u.z()) + 1);
  low = Max(low, getExtraMeshElementLowIndex());
  high = Min(high, getExtraMeshElementHighIndex());
  return MeshElementIterator(low, high);
}

/**
* This will return an iterator which will include all the meshNodes
* contained by the bounding box which also intersect the mesh.
* If a dimension of the widget is degenerate (has a thickness of 0)
* the nearest meshNode in that dimension is used.
*
* The mesh region includes extra meshNodes
*/
MeshNodeIterator
Mesh::getMeshNodeIterator(const Box& b) const
{
  // Determine if we are dealing with a 2D box.
   Point l = getLevel()->positionToIndex(b.lower());
   Point u = getLevel()->positionToIndex(b.upper());
   int low_x, low_y, low_z, high_x, high_y, high_z;
   if (l.x() != u.x()) {
     // Get the meshNodes that are included
     low_x = RoundUp(l.x());
     high_x = RoundDown(u.x()) + 1;
   } else {
     // Get the meshNodes that are nearest
     low_x = RoundDown(l.x()+0.5);
     high_x = low_x + 1;
   }
   if (l.y() != u.y()) {
     // Get the meshNodes that are included
     low_y = RoundUp(l.y());
     high_y = RoundDown(u.y()) + 1;
   } else {
     // Get the meshNodes that are nearest
     low_y = RoundDown(l.y()+0.5);
     high_y = low_y + 1;
   }
   if (l.z() != u.z()) {
     // Get the meshNodes that are included
     low_z = RoundUp(l.z());
     high_z = RoundDown(u.z()) + 1;
   } else {
     // Get the meshNodes that are nearest
     low_z = RoundDown(l.z()+0.5);
     high_z = low_z + 1;
   }
   IntVector low(low_x, low_y, low_z);
   IntVector high(high_x, high_y, high_z);
   low = Max(low, getExtraMeshNodeLowIndex());
   high = Min(high, getExtraMeshNodeHighIndex());
   return MeshNodeIterator(low, high);
}

// void Mesh::initializeBoundaryConditions()
// {
// }

