#ifndef __SMOOTH_SPHERE_PIECE_H__
#define __SMOOTH_SPHERE_PIECE_H__

#include <Core/GeometryPiece/SmoothGeomPiece.h>
#include <Core/Geometry/Point.h>

#include <cmath>
#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif

namespace Uintah {

/////////////////////////////////////////////////////////////////////////////
/*!
	
 \class SmoothSphereGeomPiece
	
 \brief Creates a smooth sphere geometry
	
 \author  Biswajit Banerjee \n
          Callaghan Innovation \n

   Creates a smooth solid/hollow sphere with/without end-caps from the 
   xml input file description.
   The input form for a solid sphere looks like this: \n
   \verbatim
   <smooth_sphere> 
     <center>[0.0,0.0,0.0]</center> 
     <outer_radius>2.0</outer_radius> 
     <num_radial_pts>20</num_radial_pts> 
   </smooth_sphere> 
   \endverbatim
   The input form for a hollow sphere looks like this: \n
   \verbatim
   <smooth_sphere> 
     <center>[0.0,0.0,0.0]</center> 
     <outer_radius>2.0</outer_radius> 
     <inner_radius>1.0</inner_radius> 
     <num_radial_pts>20</num_radial_pts> 
   </smooth_sphere> 
   \endverbatim
   If the points are to be written to an output file, use the following
   \verbatim
   <smooth_sphere> 
     <center>[0.0,0.0,0.0]</center> 
     <outer_radius>2.0</outer_radius> 
     <inner_radius>1.0</inner_radius> 
     <num_radial_pts>20</num_radial_pts> 
     <output_file>"fileName"</output_file>
   </smooth_sphere> 
   \endverbatim
	
*/
/////////////////////////////////////////////////////////////////////////////

  class SmoothSphereGeomPiece : public SmoothGeomPiece {
	 
  public:
    //////////////////////////////////////////////////////////////////////
    /*!  
      \brief Constructor that takes a ProblemSpecP argument.   
      It reads the xml input specification and builds a sphere.
    */
    //////////////////////////////////////////////////////////////////////
    SmoothSphereGeomPiece(ProblemSpecP &);
	 
    //////////////////////////////////////////////////////////////////////
    /*! Destructor */
    //////////////////////////////////////////////////////////////////////
    virtual ~SmoothSphereGeomPiece();

    static const string TYPE_NAME;
    virtual std::string getType() const { return TYPE_NAME; }

    /// Make a clone
    virtual GeometryPieceP clone() const;
	 
    //////////////////////////////////////////////////////////////////////
    /*! Determines whether a point is inside the sphere. */
    //////////////////////////////////////////////////////////////////////
    virtual bool inside(const Point &p) const;
	 
    //////////////////////////////////////////////////////////////////////
    /*! Returns the bounding box surrounding the box. */
    //////////////////////////////////////////////////////////////////////
    virtual Box getBoundingBox() const;

    //////////////////////////////////////////////////////////////////////
    /*! Creates the particles */
    //////////////////////////////////////////////////////////////////////
    virtual unsigned int createPoints();

  private:

    //////////////////////////////////////////////////////////////////////
    /*! Create points that discretize the sphere */
    //////////////////////////////////////////////////////////////////////
    int createSpherePointsSpiral();
    int createSpherePointsEqualArea();

    virtual void outputHelper( ProblemSpecP & ps ) const;

    //////////////////////////////////////////////////////////////////////
    /*! Create the point set on a unit 2-sphere with origin at (0.0, 0.0, 0.0) 
         using Leopardi's recursive algorithm */
    //////////////////////////////////////////////////////////////////////
    void createPointSetSpiral(double shell_outer_radius,
                              double shell_inner_radius,
                              double characteristic_size);
    void createPointSetPolar2D(double shell_outer_radius,
                               double shell_inner_radius,
                               double characteristic_size);

    //////////////////////////////////////////////////////////////////////
    /*! Create the point set on a unit 1-sphere with origin at (0.0, 0.0, 0.0) */
    //////////////////////////////////////////////////////////////////////
    void createPointSetPolar1D(int num_points,
                               std::vector<double>& points);
	 
    //////////////////////////////////////////////////////////////////////
    /*! Create the set of nested caps in 2D */
    //////////////////////////////////////////////////////////////////////
    void createCaps2D(int num_points,
                      std::vector<double>& cap_colatitudes,
                      std::vector<int>& int_regions);

    //////////////////////////////////////////////////////////////////////
    /*  Create the set of nested caps in 1D */
    //////////////////////////////////////////////////////////////////////
    void createCaps1D(int num_points,
                      std::vector<double>& cap_colatitudes);
    //////////////////////////////////////////////////////////////////////
    /*  Maximize the minimum distance of center points of collars */
    //////////////////////////////////////////////////////////////////////
    double circleOffset(int num_start, int num_end);

    //////////////////////////////////////////////////////////////////////
    /*  Compute the colatitide of the north polar spherical cap */
    //////////////////////////////////////////////////////////////////////
    double polarColatitude(int num_points);

    //////////////////////////////////////////////////////////////////////
    /*  Compute the numbers of collars between the polar caps */
    //////////////////////////////////////////////////////////////////////
    int numCollars(int num_points,
                   double polar_colatitude,
                   double ideal_collar_angle);

    //////////////////////////////////////////////////////////////////////
    /*  Compute ideal real number of regions in each collar */
    //////////////////////////////////////////////////////////////////////
    void idealRegionList(int num_points,
                         double polar_colatitude,
                         int num_collars,
                         std::vector<double>& real_regions);

    //////////////////////////////////////////////////////////////////////
    /*  Round to nearest int */
    //////////////////////////////////////////////////////////////////////
    void roundToNaturals(int num_points,
                         std::vector<double>& real_regions,
                         std::vector<int>& int_regions);

    //////////////////////////////////////////////////////////////////////
    /*  Compute colatitudes of spherical caps enclosing cumulative sum of regions */
    //////////////////////////////////////////////////////////////////////
    void capColatitudes(int num_points, 
                        double polar_colatitude,
                        std::vector<int>& int_regions,
                        std::vector<double>& cap_colatitudes);

    //////////////////////////////////////////////////////////////////////
    /*  Compute the spherical radius of a cap */
    //////////////////////////////////////////////////////////////////////
    inline double sRadiusOfCap(double area);

    //////////////////////////////////////////////////////////////////////
    /*  Compute angle for spherecal collars for an equal area partition */
    //////////////////////////////////////////////////////////////////////
    inline double idealCollarAngle(int num_points);

    //////////////////////////////////////////////////////////////////////
    /*  Compute area of one region of a equal area partition */
    //////////////////////////////////////////////////////////////////////
    inline double areaOfIdealRegion(int num_points);

    //////////////////////////////////////////////////////////////////////
    /*  Compute the area of a spherical collar */
    //////////////////////////////////////////////////////////////////////
    inline double areaOfCollar(double collar_start_radius, 
                               double collar_end_radius);

    //////////////////////////////////////////////////////////////////////
    /*  Compute the area of a spherical cap */
    //////////////////////////////////////////////////////////////////////
    inline double areaOfCap(double collar_radius);

    //////////////////////////////////////////////////////////////////////
    /*  Compute the surface area of a sphere  */
    //////////////////////////////////////////////////////////////////////
    inline double surfaceAreaOfSphere(double radius);

    //////////////////////////////////////////////////////////////////////
    /*  Compute the volume of a sphere */
    //////////////////////////////////////////////////////////////////////
    inline double volumeOfSphere(double radius);

    //////////////////////////////////////////////////////////////////////
    /*  Compute greatest common divisor of two ints */
    //////////////////////////////////////////////////////////////////////
    int gcd(int u, int v);

    //////////////////////////////////////////////////////////////////////
    /*  Compute complete elliptic integral of second kind
        In this case mm < 1 so imaginary modulus transformation is needed*/
    //////////////////////////////////////////////////////////////////////
    double elliptic2E(double mm);

    Point  d_center;
    double d_outerRadius;
    double d_innerRadius;
    int d_numRadial;
    string d_fileName;

  };
} // End namespace Uintah

#endif // __SMOOTH_SPHERE_PIECE_H__
