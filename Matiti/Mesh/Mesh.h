#ifndef MATITI_MESH_H
#define MATITI_MESH_H

#include <Mesh/Domain.h>
#include <Core/Domain/fixedvector.h>
#include <Mesh/MeshElementIterator.h>
#include <Mesh/MeshNodeIterator.h>
#include <Core/Domain/Variables/Iterator.h>

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Exceptions/InternalError.h>

#include   <string>
#include   <map>
#include   <iosfwd>
#include   <vector>


namespace Matiti {

  using SCIRun::Vector;
  using SCIRun::Point;
  using SCIRun::IntVector;

  class MeshNodeIterator;
  class MeshElementIterator;
   
  class Mesh {
    public:

    friend std::ostream& operator<<(std::ostream& out, const Matiti::Mesh & r);

    class Compare {
      public:
        inline bool operator()(const Mesh* p1, const Mesh* p2) const {
            return (p1 != 0 && p2 != 0) ? (p1->getID() < p2->getID()) :
              ((p2 != 0) ? true : false);
          }
      private:
    };
      
    static const int MAX_MESH_SELECT = 32; 
    typedef fixedvector<const Mesh*, MAX_MESH_SELECT> selectType;

    inline MeshElement* getMeshElementLow() const
    {
      return d_lowElement;
    }

    inline MeshElement* getMeshElementHigh() const
    {
      return d_highElement;
    }

    inline MeshNode* getMeshNodeLow() const
    {
      return d_lowNode;
    }

    inline IntVector getMeshNodeHigh() const
    {
      return d_highNode;
    } 

    inline MeshElementIterator getMeshElementIterator() const
    {
      return MeshElementIterator(getMeshElementLow(),getMeshElementHigh());
    }

    inline MeshNodeIterator getMeshNodeIterator() const
    {
      return MeshNodeIterator(getMeshNodeLow(),getMeshNodeHigh());
    }

    /**
       * Sets a pointer to the new grid
    */
    inline void setDomain(Domain* domain)
    {
      //set the domain pointer 
      d_domain=domain;
    }

    /**
       * Returns the domain coordinates of the node idx
    */
    inline Point getMeshNodePosition(/*const IntVector& idx*/) const 
    {
    }

    /**
       * Returns the domain coordinates of the cell idx
    */
    inline Point getMeshElementPosition(/*const IntVector& idx*/) const 
    {
    }
      
    /**
       * Returns the cell index of the coordinate pos
    */
    inline IntVector getMeshElementIndex(const Point& pos) const 
    {
    }

    /**
       * Returns the 8 nodes found around the point pos
    */
    void findMeshElementMeshNodes(const Point& pos, IntVector ni[8]) const;

    /**
       * Returns the 27 nodes found around the point pos
    */
    void findMeshElementMeshNodes27(const Point& pos, IntVector ni[27]) const;

    /**
       * Returns true if the point p is contained within the patch
       * excluding extra cells
    */
    inline bool containsPoint(const Point& p) const {
      return true;
    }

    /**
       * Returns the cell that contains the point pos
    */
    inline bool findMeshElement(const Point& pos, IntVector& ci) const
    {
      return true;
    }

    /**
       * sets the array cellIndex equal to the 8 cells
       * that contribute to the node nodeIndex.
    */
    static void findMeshElementsFromMeshNode( const IntVector& nodeIndex,
          IntVector cellIndex[8]);

    /**
       * sets the array nodeIndex equal to the 8 nodes that this
       * cell contributes to
    */
    static void findMeshNodesFromMeshElement( const IntVector& cellIndex,
          IntVector nodeIndex[8]);

    /**
       * Returns true if the node idx is owned by this patch
       * including extra cells
    */
    inline bool containsMeshNode(const IntVector& idx) const {
      return true;
    }

    /**
       * Returns true if the cell idx is owned by this patch
       * including extra cells
    */
    inline bool containsMeshElement(const IntVector& idx) const {
      return true;
    }

    /**
       * Returns the closest node to the point pos.  This node
       * is not guarenteed to be owned by this patch.
    */
    inline IntVector findClosestMeshNode(const Point& pos) const
    {
      return IntVector(0,0,0);
    }

    /**
       * Returns the position of the node idx in domain coordinates.
    */
    Point meshNodePosition(const IntVector& idx) const;

    /**
       * Returns the position of the cell idx in domain coordinates.
    */
    Point meshElementPosition(const IntVector& idx) const;

    /**
       * prints the patch boundary conditions to to ostream out.
    */
    void printMeshBCs(std::ostream& out) const;

    /**
       * returns a unique patch id
    */
    inline int getID() const {
      return d_id;
    }
      
    // void setBCType(FaceType face, BCType newbc);

    //  void initializeBoundaryConditions();
      
    protected:

      friend class Domain;     
      friend class MeshNodeIterator;     

      Mesh(int id=-1);
      ~Mesh();

    private:

      Domain* d_domain;

      /**
       * A unique patch id.
       */
      int d_id;

      /** 
       * A pointer to the real patch 
       * This pointer is null for non-virtual patches.
       * Virtual patches exist for wrap-around on periodic boundary 
       * conditions.
       */
      const Mesh* d_mesh;

      Mesh(const Mesh&);
      Mesh(const Mesh* mesh);
      Mesh& operator=(const Mesh&);
     
   }; // end class Mesh

} // End namespace Matiti

#endif
