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

#include <Core/Grid/Variables/NeighborConnectivity.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Util/Endian.h>
#include <Core/Util/FancyAssert.h>
#include <Core/Malloc/Allocator.h>


namespace Uintah {

  std::ostream& operator<<(std::ostream &out, 
                           const Uintah::NeighborConnectivity& conn) 
  {
    for (int ii = 0; ii < 216; ii++) {
      out << conn.d_connected[ii] << " ";
    }
    return out;
  }


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
