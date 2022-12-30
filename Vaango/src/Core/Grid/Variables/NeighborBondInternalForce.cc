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

#include <Core/Grid/Variables/NeighborBondInternalForce.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Util/Endian.h>
#include <Core/Util/FancyAssert.h>
#include <Core/Malloc/Allocator.h>

using namespace Uintah;

const std::string& 
NeighborBondInternalForce::get_h_file_path()
{
  static const std::string path(Uintah::FETypeDescription::cc_to_h(__FILE__));
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
namespace Uintah {

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

  const FETypeDescription* 
  get_fetype_description(Uintah::NeighborBondInternalForce*)
  {
    static FETypeDescription* td = 0;
    if (!td) {
      td = scinew FETypeDescription("NeighborBondInternalForce", 
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

} // namespace Uintah

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
      td = scinew TypeDescription(TypeDescription::Type::NeighborBondInternalForce, "NeighborBondInternalForce", true,
                                  &makeMPI_NeighborBondInternalForce);
    }
    return td;
  }

} // End namespace Uintah
