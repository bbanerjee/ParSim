#ifndef __SPHEREGEOMETRYPIECE_H__
#define __SPHEREGEOMETRYPIECE_H__

#include <GeometryPiece/GeometryPiece.h>
#include <GeometryMath/Box3D.h>
#include <GeometryMath/Point3D.h>
#include <GeometryMath/Vector3D.h>
#include <GeometryMath/IntVector3D.h>
#include <MPMPatch.h>
#include <MPMDataTypes.h>
//#include <MPMData.h>

#include <Core/ProblemSpec/ProblemSpecP.h>


namespace BrMPM {
     
     class SphereGeometryPiece : public GeometryPiece
     {
      public:
       
       SphereGeometryPiece();

       virtual ~SphereGeometryPiece();
       
       void
       initialise(MPMPatch& Patch);

       SphereGeometryPiece(Uintah::ProblemSpecP& ps,
                           Point3DParticleData& setPointsOfSphere);

       void outerBox(Box3D& Box);

       bool inside (const Point3D& pt) const;

       std::string name();

      protected:

       void createParticles(Point3DParticleData& setPoints);

      private:

       Vector3D d_ghost;
       IntVector3D d_num_particles_per_cell;

       double d_radius;
       
       MPMPatch d_patch;

       Box3D d_outerbox;

       Point3D d_center;
       Point3D d_outerlower;
       Point3D d_outerupper;
       

      // Vector3D d_node_counts;
       Vector3D d_num_grids;
       Vector3D d_cellsize;
       
    
       //int d_t_initial;
       //int d_t_final;





       //std::string d_name;

       //double xcoord, ycoord, zcoord;



      }; //end of class

} //end of BrMPM namespace
       
          


































#endif  //ifndef __SPHEREGEOMETRYPIECE_H__
