#ifndef ELLIP3D_BOUNDARY_FILE_READER_H_H
#define ELLIP3D_BOUNDARY_FILE_READER_H_H

#include <Boundary/BoundaryContainers.h>
#include <Core/Geometry/Box.h>
#include <Core/Geometry/OrientedBox.h>
#include <Core/Types/RealTypes.h>

namespace dem {

class BoundaryFileReader
{

public:
  BoundaryFileReader() = default;
  ~BoundaryFileReader() = default;

  void read(const std::string& inputFileName, Box& domain, Box& patchBox,
            BoundaryPArray& boundaries) const;
  bool readXML(const std::string& inputFileName, Box& domain, Box& patchBox,
               BoundaryPArray& boundaries) const;
  bool readJSON(const std::string& inputFileName, Box& domain, Box& patchBox,
                BoundaryPArray& boundaries) const;

  void readVTK(const std::string& filename, OrientedBox& domain) const;

private:

  BoundaryFileReader(BoundaryFileReader const&) = delete; // don't implement
  void operator=(BoundaryFileReader const&) = delete; // don't implement
};
}
#endif
