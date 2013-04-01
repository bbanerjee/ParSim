#ifndef MATITI_MESH_H
#define MATITI_MESH_H

#include <Vaango/Core/Geometry/Vector.h>
#include <Mesh/MeshNode.h>
#include <Mesh/MeshBond.h>
#include <vector>

namespace Matiti {

  class Mesh
  {
    public:
      Mesh();
      ~Mesh();

    private:

      std::vector<MeshNode> d_nodes;
      std::vector<MeshBond> d_bonds;

      std::vector<SCIRun::Vector> d_displacement;
      std::vector<SCIRun::Vector> d_velocity;
      std::vector<SCIRun::Vector> d_acceleration;
      std::vector<SCIRun::Vector> d_force;

      // prevent copying
      Mesh(const Mesh& mesh);
      Mesh& operator=(const Mesh& mesh);

  };
} // End namespace Matiti

#endif // MATITI_MESH_H
