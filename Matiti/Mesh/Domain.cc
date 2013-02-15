#include <Mesh/Domain.h>
#include <Mesh/Mesh.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Math/MatitiMiscMath.h>
#include <Core/Math/Primes.h>
#include <Core/Math/MiscMath.h>
#include <iostream>
#include <sci_values.h>

using namespace std;
using namespace Matiti;

Domain::Domain()
{
  all_meshes = 0;
}

Domain::~Domain()
{
  for (meshIterator iter=d_meshes.begin(); iter != d_meshes.edn(); iter++) {
    delete *iter;
  }
  if (all_meshes && all_meshes->removeReference()) delete all_meshes;
}

void 
Domain::problemSetup(const Uintah::ProblemSpecP& params)
{
   Uintah::ProblemSpecP domain_ps = params->findBlock("Domain");
   if(!domain_ps) return;

   // find upper/lower corner
   Uintah::ProblemSpecP box_ps = domain_ps->findBlock("Box");
   box_ps->require("lower", d_lower);
   box_ps->require("upper", d_upper);

} // end problemSetup()

Uintah::Point Domain::getUpper() const
{
  return d_upper;
}

Uintah::Point Domain::getLower() const
{
  return d_lower;
}

bool Uintah::containsPoint(const Uintah::Point& pt) const
{
  Uintah::Vector vlo(pt-d_lower);
  Uintah::Vector vhi(pt-d_upper);
  return vlo.x() >= 0.0 && vlo.y() >= 0.0 && vlo.z() >= 0.0  
         &&  vhi.x() <= 0.0 && vhi.y() <= 0.0 && vhi.z() <= 0.0 ; 
}

Domain::const_meshIterator Domain::meshesBegin() const
{
  return d_meshes.begin();
}

Domain::const_meshIterator Domain::meshesEnd() const
{
  return d_meshes.end();
}

Domain::meshIterator Domain::meshesBegin()
{
  return d_meshes.begin();
}

Domain::meshIterator Domain::meshesEnd()
{
  return d_meshes.end();
}

Mesh* Domain::addMesh() 
{
  Mesh* mesh = new Mesh(this);
  d_meshes.push_back(mesh);
  return mesh;
}

int Domain::numMeshes() 
{
  return (int) d_meshes.size();
}

const Mesh* Domain::getMeshFromPoint(const Uintah::Point& ) const
{
  std::cerr << "getMeshFromPoint not implemeneted yet" << std::endl;
}

const MeshSet* Domain::allMeshes() const
{
  return all_meshes;
}

void Domain::finalizeDomain()
{
  vector<const Mesh*> tmp_meshes(d_meshes.size());
  for (int ii=0; ii < (int) d_meshes.size(); ii++) {
    tmp_meshes[ii] = d_meshes[ii];
  }
  all_meshes = new MeshSet();
  all_meshes->addReference();
  all_meshes->addAll(tmp_meshes);

  for (meshIterator iter=d_meshes.begin(); iter != d_meshes.end(); iter++) {
    (*iter)->finalizeMesh();
  }
}

void Domain::assignBCS(const Uintah::ProblemSpecP &domain_ps)
{
  Uintah::ProblemSpecP bc_ps = domain_ps->findBlock("BoundaryConditions");
  if (bc_ps == 0) return;

  //BoundaryConditionReader reader;
  //reader.read(bc_ps);

  //for (Domain::FaceType face_side = Domain::startFace;
  //     face_side <= Domain::endFace; face_side=Domain::nextFace(face_side)) {
  //  if (getBCType(face_side) == Domain::None) {
  //    setArrayBCValues(face_side, &(reader.d_BCReaderData[face_side]));
  //  }
  //}
}

//void
//Patch::setArrayBCValues(Patch::FaceType face, BCDataArray* bc)
//{
//  // At this point need to set up the iterators for each BCData type:
//  // Side, Rectangle, Circle, Difference, and Union.
//
//  bc->determineIteratorLimits(face,this);
//  (*d_arrayBCS)[face] = bc->clone();
//}


namespace Matiti
{
  ostream& operator<<(ostream& out, const Domain& domain)
  {
    out.setf(ios::floatfield);
    out.precision(6);
    out << "Domain: " << domain.getLower() << "," << domain.getUpper() << endl;
    out << " has " << numMeshes() << " mesh(es)" << endl;
    for ( Level::meshIterator meshIter = meshesBegin(); meshIter < meshesEnd(); meshIter++ ) {
      const Mesh* mesh = *meshIter;
      out <<"    "<< *mesh << endl;
    }
    return out;
  }
}

