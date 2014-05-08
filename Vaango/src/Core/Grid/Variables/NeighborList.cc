#include <Core/Grid/Variables/NeighborList.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Util/Endian.h>
#include <Core/Util/FancyAssert.h>
#include <Core/Malloc/Allocator.h>

using namespace Uintah;

const std::string& 
NeighborList::get_h_file_path()
{
  static const std::string path(SCIRun::TypeDescription::cc_to_h(__FILE__));
  return path;
}

// Added for compatibility with core types
namespace SCIRun {

  void 
  swapbytes(Uintah::NeighborList& family)
  {
    Uintah::ParticleID* ptr = (Uintah::ParticleID*) (&family);
    SWAP_8(*ptr);
    for (int ii = 1; ii < 216; ii++) {
      SWAP_8(*++ptr);
    }
  }

  template<> const std::string 
  find_type_name(Uintah::NeighborList*)
  {
    static const std::string name = "NeighborList";
    return name;
  }

  const TypeDescription* 
  get_type_description(Uintah::NeighborList*)
  {
    static TypeDescription* td = 0;
    if (!td) {
      td = scinew TypeDescription("NeighborList", 
                                  Uintah::NeighborList::get_h_file_path(),
                                  "Uintah");
    }
    return td;
  }

  void 
  Pio(Piostream& stream, Uintah::NeighborList& family)
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
  MPI_Datatype makeMPI_NeighborList()
  {
    ASSERTEQ(sizeof(NeighborList), sizeof(ParticleID)*216);

    MPI_Datatype mpitype;
    MPI_Type_vector(1, 216, 216, MPI_LONG_LONG_INT, &mpitype);
    MPI_Type_commit(&mpitype);

    return mpitype;
  }

  const TypeDescription* fun_getTypeDescription(NeighborList*)
  {
    static TypeDescription* td = 0;
    if(!td){
      td = scinew TypeDescription(TypeDescription::NeighborList, "NeighborList", true,
                                  &makeMPI_NeighborList);
    }
    return td;
  }

} // End namespace Uintah
