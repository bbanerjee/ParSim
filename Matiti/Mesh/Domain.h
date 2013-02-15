#ifndef MATITI_DOMAIN_H
#define MATITI_DOMAIN_H

#include <Mesh/DomainP.h>
#include <Common/Handle.h>
#include <Common/RefCounted.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <vector>
#include <list>

namespace Matiti {

  class Domain : public RefCounted {
  public:
    typedef std::vector<Mesh*>::iterator meshIterator;
    typedef std::vector<Mesh*>::const_iterator const_meshIterator;

    friend std::ostream& operator<<(std::ostream& out, const Matiti::Domain& mesh);

    enum BCType {
      None=0,
      Symmetry=1
    };

    enum FaceType {
      xminus = 0,
      xplus = 1,
      yminus = 2,
      yplus = 3,
      zminus = 4,
      zplus = 5,
      startFace = xminus,
      endFace = zplus,
      numFaces,
      invalidFace
    };

    Domain();
    virtual ~Domain();

    // Problem setup functions called from simulation controller
    void problemSetup(const ProblemSpecP& params);
    Uintah::Point getUpper() const;
    Uintah::Point getLower() const;
    bool containsPoint(const Uintah::Point& ) const;

    const_meshIterator meshesBegin() const;
    const_meshIterator meshesEnd() const;

    meshIterator meshesBegin();
    meshIterator meshesEnd();

    Mesh* addMesh();
    int numMeshes();

    const Mesh* getMesh(int index) const {return d_meshes[index];}
    const Mesh* getMeshFromPoint(const Uintah::Point&) const;
    const MeshSet* allMeshes() const;

    void finalizeDomain();

    //Assigns the boundary conditions to the domain
    void assignBCS( const ProblemSpecP &mesh_ps);

    static inline FaceType nextFace(FaceType face) {
      return (FaceType)((int)face+1);
    }

    inline BCType getBCType(FaceType face) const {
      switch(face)
      {
        case xminus:
          return static_cast<BCType>(d_domainState.xminus);
        case yminus:
          return static_cast<BCType>(d_domainState.yminus);
        case zminus:
          return static_cast<BCType>(d_domainState.zminus);
        case xplus:
          return static_cast<BCType>(d_domainState.xplus);
        case yplus:
          return static_cast<BCType>(d_domainState.yplus);
        case zplus:
          return static_cast<BCType>(d_domainState.zplus);
        default:
          throw SCIRun::InternalError("Invalid FaceType Specified", __FILE__, __LINE__);
          return None;
      }
    }

    //void setArrayBCValues(FaceType face, BCDataArray* bc);

  private:
    
    Domain(const Domain&);
    Domain& operator=(const Domain&);

    MeshSet* all_meshes;
    std::vector<Mesh*> d_meshes;
    Uintah::Point d_upper, d_lower;
    std::vector<BCDataArray*>* d_arrayBCs;

    struct DomainState
    {
      //The boundary conditions for each face
      unsigned int xminus : 2;
      unsigned int xplus : 2;
      unsigned int yminus : 2;
      unsigned int yplus : 2;
      unsigned int zminus : 2;
      unsigned int zplus : 2;
    };
    DomainState d_domainState;
  };

} // End namespace Matiti

#endif
