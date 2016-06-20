/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <Core/Grid/Variables/NeighborList.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Util/Endian.h>
#include <Core/Util/FancyAssert.h>
#include <Core/Malloc/Allocator.h>

using namespace Uintah;

const std::string& 
NeighborList::get_h_file_path()
{
  static const std::string path(Uintah::TypeDescription::cc_to_h(__FILE__));
  return path;
}

namespace Uintah {

  std::ostream& operator<<(std::ostream &out, 
                           const Uintah::NeighborList& family) 
  {
    for (int ii = 0; ii < 216; ii++) {
      out << family.d_family[ii] << " ";
    }
    return out;
  }
}

// Added for compatibility with core types
namespace Uintah {

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

} // namespace Uintah

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
