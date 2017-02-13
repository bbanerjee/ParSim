#ifndef ELLIP3D_BOUNDARY_READER_H_H
#define ELLIP3D_BOUNDARY_READER_H_H

#include <Core/Types/realtypes.h>
#include <Core/Geometry/Box.h>

namespace dem {

class BoundaryReader
{

public:

  BoundaryReader() = default;
  ~BoundaryReader() = default;

  void read(const char* input,
            Box& container,
            Box& grid,
            BoundaryPArray& boundaries) const;
  bool readXML(const std::string& inputFileName,
               Box& container,
               Box& grid,
               BoundaryPArray& boundaries) const;

private:

  enum BoundaryTypes {
    PLANE,
    CYLINDER,
    NONE
  };

  BoundaryTypes getEnum(const std::string& str) const {
    if (str == "plane") return BoundaryTypes::PLANE;
    else if (str == "cylinder") return BoundaryTypes::CYLINDER;
    else return BoundaryTypes::NONE;
  }

  // make sure these two are unaccessable to avoid copies of singelton
  BoundaryReader(BoundaryReader const&) = delete;      // don't implement
  void operator=(BoundaryReader const&) = delete; // don't implement

};
}
#endif
