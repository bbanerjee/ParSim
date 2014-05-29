#include <Core/Grid/Variables/NeighborBondInternalForce.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Util/Endian.h>
#include <Core/Util/FancyAssert.h>
#include <Core/Malloc/Allocator.h>

using namespace Uintah;

const std::string& 
NeighborBondInternalForce::get_h_file_path()
{
  static const std::string path(SCIRun::TypeDescription::cc_to_h(__FILE__));
  return path;
}

namespace Uintah {

  std::ostream& operator<<(std::ostream &out, 
                           const Uintah::NeighborBondInternalForce& force) 
  {
    out.setf(std::ios::floatfield);
    out.precision(3);
    for (int ii = 0; ii < 216; ii++) {
      out << force.d_bondInternalForce[ii] << " ";
    }
    return out;
  }
}

// Added for compatibility with core types
namespace SCIRun {

  void 
  swapbytes(Uintah::NeighborBondInternalForce& bondInternalForce)
  {
    for (int ii = 1; ii < 216; ii++) {
      swapbytes(bondInternalForce[ii]);
    }
  }

  template<> const std::string 
  find_type_name(Uintah::NeighborBondInternalForce*)
  {
    static const std::string name = "NeighborBondInternalForce";
    return name;
  }

  const TypeDescription* 
  get_type_description(Uintah::NeighborBondInternalForce*)
  {
    static TypeDescription* td = 0;
    if (!td) {
      td = scinew TypeDescription("NeighborBondInternalForce", 
                                  Uintah::NeighborBondInternalForce::get_h_file_path(),
                                  "Uintah");
    }
    return td;
  }

  void 
  Pio(Piostream& stream, Uintah::NeighborBondInternalForce& bondInternalForce)
  {
    stream.begin_cheap_delim();
    for (int ii = 0; ii < 216; ii++) {
      Pio(stream, bondInternalForce[ii]);
    }
    stream.end_cheap_delim();
  }

} // namespace SCIRun

namespace Uintah {
  //* TODO: Serialize **/
  MPI_Datatype makeMPI_NeighborBondInternalForce()
  {
    ASSERTEQ(sizeof(NeighborBondInternalForce), sizeof(double)*(3*216));

    MPI_Datatype mpitype;
    MPI_Type_vector(1, 3*216, 3*216, MPI_DOUBLE, &mpitype);
    MPI_Type_commit(&mpitype);

    return mpitype;
  }

  const TypeDescription* fun_getTypeDescription(NeighborBondInternalForce*)
  {
    static TypeDescription* td = 0;
    if(!td){
      td = scinew TypeDescription(TypeDescription::NeighborBondInternalForce, "NeighborBondInternalForce", true,
                                  &makeMPI_NeighborBondInternalForce);
    }
    return td;
  }

} // End namespace Uintah
