#ifndef ELLIP3D_BOUNDARY_READER_H_H
#define ELLIP3D_BOUNDARY_READER_H_H

#include <Boundary/BoundaryContainers.h>
#include <Core/Geometry/Box.h>
#include <Core/Types/RealTypes.h>

namespace dem {

class BoundaryReader
{

public:
  BoundaryReader() = default;
  ~BoundaryReader() = default;

  void read(const std::string& inputFileName, Box& container, Box& patchBox,
            BoundaryPArray& boundaries) const;
  bool readXML(const std::string& inputFileName, Box& container, Box& patchBox,
               BoundaryPArray& boundaries) const;
  bool readJSON(const std::string& inputFileName, Box& container, Box& patchBox,
                BoundaryPArray& boundaries) const;

private:

  BoundaryReader(BoundaryReader const&) = delete; // don't implement
  void operator=(BoundaryReader const&) = delete; // don't implement
};
}
#endif
