#ifndef __MATITI_GEOMETRY_PIECE_H__
#define __MATITI_GEOMETRY_PIECE_H__

#include <Geometry/Box3D.h>
#include <Geometry/Point3D.h>

#include <string>

namespace BrMPM 
{

  class GeometryPiece
  {
  public:

    GeometryPiece();
    virtual ~GeometryPiece();

    virtual Box3D boundingBox() const = 0;  

    virtual bool inside (const Point3D& pt) const = 0;

    virtual std::string name() const = 0;

    void name(const std::string& name) {d_hasName = true; d_name = name;}

  protected:

    bool d_hasName;
    std::string d_name;

  }; // end class

} // end namespace

#endif
