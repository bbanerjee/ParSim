#include <Core/Grid/Variables/NeighborBondEnergy.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Util/Endian.h>
#include <Core/Util/FancyAssert.h>
#include <Core/Malloc/Allocator.h>

using namespace Uintah;

const std::string& 
NeighborBondEnergy::get_h_file_path()
{
  static const std::string path(SCIRun::TypeDescription::cc_to_h(__FILE__));
  return path;
}

// Added for compatibility with core types
namespace SCIRun {

  void 
  swapbytes(Uintah::NeighborBondEnergy& family)
  {
    double* ptr = (double*) (&family);
    SWAP_8(*ptr);
    for (int ii = 1; ii < 216; ii++) {
      SWAP_8(*++ptr);
    }
  }

  template<> const std::string 
  find_type_name(Uintah::NeighborBondEnergy*)
  {
    static const std::string name = "NeighborBondEnergy";
    return name;
  }

  const TypeDescription* 
  get_type_description(Uintah::NeighborBondEnergy*)
  {
    static TypeDescription* td = 0;
    if (!td) {
      td = scinew TypeDescription("NeighborBondEnergy", 
                                  Uintah::NeighborBondEnergy::get_h_file_path(),
                                  "Uintah");
    }
    return td;
  }

  void 
  Pio(Piostream& stream, Uintah::NeighborBondEnergy& family)
  {
    stream.begin_cheap_delim();
    for (int ii = 0; ii < 216; ii++) {
      Pio(stream, family[ii]);
    }
    stream.end_cheap_delim();
  }

} // namespace SCIRun

namespace Uintah {
  //* TODO: Serialize **/
  MPI_Datatype makeMPI_NeighborBondEnergy()
  {
    ASSERTEQ(sizeof(NeighborBondEnergy), sizeof(double)*216);

    MPI_Datatype mpitype;
    MPI_Type_vector(1, 216, 216, MPI_DOUBLE, &mpitype);
    MPI_Type_commit(&mpitype);

    return mpitype;
  }

  const TypeDescription* fun_getTypeDescription(NeighborBondEnergy*)
  {
    static TypeDescription* td = 0;
    if(!td){
      td = scinew TypeDescription(TypeDescription::NeighborBondEnergy, "NeighborBondEnergy", true,
                                  &makeMPI_NeighborBondEnergy);
    }
    return td;
  }

} // End namespace Uintah
