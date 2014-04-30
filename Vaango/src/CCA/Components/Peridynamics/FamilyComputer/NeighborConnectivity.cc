#include <CCA/Components/Peridynamics/NeighborConnectivity.h>
#include <Core/Util/Endian.h>
#include <Core/Util/TypeDescription.h>

const std::string& 
NeighborConnectivity::get_h_file_path()
{
  static const std::string path(SCIRun::TypeDescription::cc_to_h(__FILE__));
  return path;
}

// Added for compatibility with core types
namespace SCIRun {

  void 
  swapbytes(Vaango::NeighborConnectivity& )
  {
    // Nothing to be done here
  }

  template<> const std::string 
  find_type_name(Vaango::NeighborConnectivity*)
  {
    static const std::string name = "NeighborConnectivity";
    return name;
  }

  const TypeDescription* 
  get_type_description(Vaango::NeighborConnectivity*)
  {
    static TypeDescription* td = 0;
    if (!td) {
      td = scinew TypeDescription("NeighborConnectivity", 
                                  Vaango::NeighborConnectivity::get_h_file_path(),
                                  "Vaango");
    }
    return td;
  }

  void 
  Pio(Piostream& stream, Vaango::NeighborConnectivity& broken)
  {
    stream.begin_cheap_delim();
    for (int ii = 0; ii < 216; ii++) {
      Pio(stream, broken[ii]);
    }
    stream.end_cheap_delim();
  }

} // namespace SCIRun


#endif
