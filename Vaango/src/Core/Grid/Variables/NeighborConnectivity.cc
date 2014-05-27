#include <Core/Grid/Variables/NeighborConnectivity.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Util/Endian.h>
#include <Core/Util/FancyAssert.h>
#include <Core/Malloc/Allocator.h>

namespace Uintah {

  const std::string& 
  NeighborConnectivity::get_h_file_path()
  {
    static const std::string path(SCIRun::TypeDescription::cc_to_h(__FILE__));
    return path;
  }
}

namespace Uintah {

  std::ostream& operator<<(std::ostream &out, 
                           const Uintah::NeighborConnectivity& conn) 
  {
    for (int ii = 0; ii < 216; ii++) {
      out << conn.d_connected[ii] << " ";
    }
    return out;
  }
}

// Added for compatibility with core types
namespace SCIRun {

  void 
  swapbytes(Uintah::NeighborConnectivity& )
  {
    // Nothing to be done here
  }

  template<> const std::string 
  find_type_name(Uintah::NeighborConnectivity*)
  {
    static const std::string name = "NeighborConnectivity";
    return name;
  }

  const TypeDescription* 
  get_type_description(Uintah::NeighborConnectivity*)
  {
    static TypeDescription* td = 0;
    if (!td) {
      td = scinew TypeDescription("NeighborConnectivity", 
                                  Uintah::NeighborConnectivity::get_h_file_path(),
                                  "Uintah");
    }
    return td;
  }

  void 
  Pio(Piostream& stream, Uintah::NeighborConnectivity& broken)
  {
    stream.begin_cheap_delim();
    for (int ii = 0; ii < 216; ii++) {
      Pio(stream, broken[ii]);
    }
    stream.end_cheap_delim();
  }

} // namespace SCIRun


namespace Uintah {
  //* TODO: Serialize **/
  MPI_Datatype makeMPI_NeighborConnectivity()
  {
    ASSERTEQ(sizeof(NeighborConnectivity), sizeof(bool)*216);

    MPI_Datatype mpitype;
    MPI_Type_vector(1, 216, 216, MPI_UB, &mpitype);
    MPI_Type_commit(&mpitype);

    return mpitype;
  }

  const TypeDescription* fun_getTypeDescription(NeighborConnectivity*)
  {
    static TypeDescription* td = 0;
    if(!td){
      td = scinew TypeDescription(TypeDescription::NeighborConnectivity, "NeighborConnectivity", true,
                                  &makeMPI_NeighborConnectivity);
    }
    return td;
  }

} // End namespace Uintah
