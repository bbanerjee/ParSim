/*
 * The MIT License
 *
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#include <Core/Disclosure/TypeDescription.h>
#include <Core/Grid/Variables/MPMIntVarTypes.h>
#include <Core/Util/Assert.h>
#include <Core/Util/Endian.h>
#include <Core/Util/FancyAssert.h>

namespace Uintah {

void
swapbytes(MetalIntVar& mp)
{
  double* p = (double*)(&mp);
  SWAP_8(*p);
  SWAP_8(*++p);
}

template<>
const std::string
find_type_name([[maybe_unused]] MetalIntVar* mp)
{
  static const std::string name = "MetalIntVar";
  return name;
}

void
swapbytes(DStressDMetalIntVar& mp)
{
  double* p = (double*)(&mp);
  SWAP_8(*p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
}

template<>
const std::string
find_type_name([[maybe_unused]] DStressDMetalIntVar* mp)
{
  static const std::string name = "DStressDMetalIntVar";
  return name;
}

void
swapbytes(ArenaIntVar& mp)
{
  double* p = (double*)(&mp);
  SWAP_8(*p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  for (int ii = 0; ii < 9; ++ii) {
    SWAP_8(*++p);
  }
}

template<>
const std::string
find_type_name([[maybe_unused]] ArenaIntVar* mp)
{
  static const std::string name = "ArenaIntVar";
  return name;
}

void
swapbytes(BorjaIntVar& mp)
{
  double* p = (double*)(&mp);
  SWAP_8(*p);
}

template<>
const std::string
find_type_name([[maybe_unused]] BorjaIntVar* mp)
{
  static const std::string name = "BorjaIntVar";
  return name;
}

void
swapbytes(SoilBrannonIntVar& mp)
{
  double* p = (double*)(&mp);
  SWAP_8(*p);
}

template<>
const std::string
find_type_name([[maybe_unused]] SoilBrannonIntVar* mp)
{
  static const std::string name = "SoilBrannonIntVar";
  return name;
}

void
swapbytes(TabularCapIntVar& mp)
{
  double* p = (double*)(&mp);
  SWAP_8(*p);
}

template<>
const std::string
find_type_name([[maybe_unused]] TabularCapIntVar* mp)
{
  static const std::string name = "TabularCapIntVar";
  return name;
}

} // end namespace Uintah

namespace Uintah {

static MPI_Datatype
makeMPI_MetalIntVar()
{
  int num_var = 2;
  ASSERTEQ(sizeof(MetalIntVar), sizeof(double) * num_var);
  MPI_Datatype mpitype;
  MPI_Type_vector(1, num_var, num_var, MPI_DOUBLE, &mpitype);
  MPI_Type_commit(&mpitype);
  return mpitype;
}

const TypeDescription*
fun_getTypeDescription(MetalIntVar*)
{
  static TypeDescription* td;
  if (!td) {
    td = scinew TypeDescription(TypeDescription::Type::Other,
                                "MetalIntVar",
                                true,
                                &makeMPI_MetalIntVar);
  }
  return td;
}

static MPI_Datatype
makeMPI_DStressDMetalIntVar()
{
  int num_var = 18;
  ASSERTEQ(sizeof(DStressDMetalIntVar), sizeof(double) * num_var);
  MPI_Datatype mpitype;
  MPI_Type_vector(1, num_var, num_var, MPI_DOUBLE, &mpitype);
  MPI_Type_commit(&mpitype);
  return mpitype;
}

const TypeDescription*
fun_getTypeDescription(DStressDMetalIntVar*)
{
  static TypeDescription* td;
  if (!td) {
    td = scinew TypeDescription(TypeDescription::Type::Other,
                                "DStressDMetalIntVar",
                                true,
                                &makeMPI_DStressDMetalIntVar);
  }
  return td;
}

static MPI_Datatype
makeMPI_ArenaIntVar()
{
  int num_var = 13;
  ASSERTEQ(sizeof(ArenaIntVar), sizeof(double) * num_var);
  MPI_Datatype mpitype;
  MPI_Type_vector(1, num_var, num_var, MPI_DOUBLE, &mpitype);
  MPI_Type_commit(&mpitype);
  return mpitype;
}

const TypeDescription*
fun_getTypeDescription(ArenaIntVar*)
{
  static TypeDescription* td;
  if (!td) {
    td = scinew TypeDescription(TypeDescription::Type::Other,
                                "ArenaIntVar",
                                true,
                                &makeMPI_ArenaIntVar);
  }
  return td;
}

static MPI_Datatype
makeMPI_BorjaIntVar()
{
  int num_var = 1;
  ASSERTEQ(sizeof(BorjaIntVar), sizeof(double) * num_var);
  MPI_Datatype mpitype;
  MPI_Type_vector(1, num_var, num_var, MPI_DOUBLE, &mpitype);
  MPI_Type_commit(&mpitype);
  return mpitype;
}

const TypeDescription*
fun_getTypeDescription(BorjaIntVar*)
{
  static TypeDescription* td;
  if (!td) {
    td = scinew TypeDescription(TypeDescription::Type::Other,
                                "BorjaIntVar",
                                true,
                                &makeMPI_BorjaIntVar);
  }
  return td;
}

static MPI_Datatype
makeMPI_SoilBrannonIntVar()
{
  int num_var = 1;
  ASSERTEQ(sizeof(SoilBrannonIntVar), sizeof(double) * num_var);
  MPI_Datatype mpitype;
  MPI_Type_vector(1, num_var, num_var, MPI_DOUBLE, &mpitype);
  MPI_Type_commit(&mpitype);
  return mpitype;
}

const TypeDescription*
fun_getTypeDescription(SoilBrannonIntVar*)
{
  static TypeDescription* td;
  if (!td) {
    td = scinew TypeDescription(TypeDescription::Type::Other,
                                "SoilBrannonIntVar",
                                true,
                                &makeMPI_SoilBrannonIntVar);
  }
  return td;
}

static MPI_Datatype
makeMPI_TabularCapIntVar()
{
  int num_var = 1;
  ASSERTEQ(sizeof(TabularCapIntVar), sizeof(double) * num_var);
  MPI_Datatype mpitype;
  MPI_Type_vector(1, num_var, num_var, MPI_DOUBLE, &mpitype);
  MPI_Type_commit(&mpitype);
  return mpitype;
}

const TypeDescription*
fun_getTypeDescription(TabularCapIntVar*)
{
  static TypeDescription* td;
  if (!td) {
    td = scinew TypeDescription(TypeDescription::Type::Other,
                                "TabularCapIntVar",
                                true,
                                &makeMPI_TabularCapIntVar);
  }
  return td;
}

static MPI_Datatype
makeMPI_ViscoScramStateData()
{
  ASSERTEQ(sizeof(ViscoScramStateData), sizeof(double) * 45);
  MPI_Datatype mpitype;
  MPI_Type_vector(1, 45, 45, MPI_DOUBLE, &mpitype);
  MPI_Type_commit(&mpitype);
  return mpitype;
}

const Uintah::TypeDescription*
fun_getTypeDescription(ViscoScramStateData*)
{
  static Uintah::TypeDescription* td = nullptr;
  if (!td) {
    td = scinew Uintah::TypeDescription(TypeDescription::Type::Other,
                                        "ViscoScramStateData",
                                        true,
                                        &makeMPI_ViscoScramStateData);
  }
  return td;
}

} // end namespace Uintah
