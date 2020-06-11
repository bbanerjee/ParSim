/*
 * The MIT License
 *
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#include <Core/Grid/Variables/MPMIntVarTypes.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Util/Endian.h>


namespace Uintah {

void swapbytes(MetalIntVar& mp) 
{
  double *p = (double *)(&mp);
  SWAP_8(*p);
  SWAP_8(*++p);
}

template <>
const std::string 
find_type_name(MetalIntVar* mp)
{
  static const std::string name="MetalIntVar";
  return name;
}

void swapbytes(ArenaIntVar& mp) 
{
  double *p = (double *)(&mp);
  SWAP_8(*p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  SWAP_8(*++p);
  for (int ii = 0; ii < 9; ++ii) {
    SWAP_8(*++p);
  }
}

template <>
const std::string 
find_type_name(ArenaIntVar* mp)
{
  static const std::string name="ArenaIntVar";
  return name;
}

void swapbytes(BorjaIntVar& mp) 
{
  double *p = (double *)(&mp);
  SWAP_8(*p);
}

template <>
const std::string 
find_type_name(BorjaIntVar* mp)
{
  static const std::string name="BorjaIntVar";
  return name;
}

void swapbytes(SoilBrannonIntVar& mp) 
{
  double *p = (double *)(&mp);
  SWAP_8(*p);
}

template <>
const std::string 
find_type_name(SoilBrannonIntVar* mp)
{
  static const std::string name="SoilBrannonIntVar";
  return name;
}

void swapbytes(TabularCapIntVar& mp) 
{
  double *p = (double *)(&mp);
  SWAP_8(*p);
}

template <>
const std::string 
find_type_name(TabularCapIntVar* mp)
{
  static const std::string name="TabularCapIntVar";
  return name;
}

} // end namespace Uintah

namespace Uintah {

static MPI_Datatype 
makeMPI_MetalIntVar()
{
  int num_var = 2;
  ASSERTEQ(sizeof(MetalIntVar), sizeof(double)*num_var);
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
    td = scinew TypeDescription(TypeDescription::Other,
                                "MetalIntVar", true,
                                &makeMPI_MetalIntVar);
  }
  return td;
}

static MPI_Datatype 
makeMPI_ArenaIntVar()
{
  int num_var = 13;
  ASSERTEQ(sizeof(ArenaIntVar), sizeof(double)*num_var);
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
    td = scinew TypeDescription(TypeDescription::Other,
                                "ArenaIntVar", true,
                                &makeMPI_ArenaIntVar);
  }
  return td;
}

static MPI_Datatype 
makeMPI_BorjaIntVar()
{
  int num_var = 1;
  ASSERTEQ(sizeof(BorjaIntVar), sizeof(double)*num_var);
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
    td = scinew TypeDescription(TypeDescription::Other,
                                "BorjaIntVar", true,
                                &makeMPI_BorjaIntVar);
  }
  return td;
}

static MPI_Datatype 
makeMPI_SoilBrannonIntVar()
{
  int num_var = 1;
  ASSERTEQ(sizeof(SoilBrannonIntVar), sizeof(double)*num_var);
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
    td = scinew TypeDescription(TypeDescription::Other,
                                "SoilBrannonIntVar", true,
                                &makeMPI_SoilBrannonIntVar);
  }
  return td;
}

static MPI_Datatype 
makeMPI_TabularCapIntVar()
{
  int num_var = 1;
  ASSERTEQ(sizeof(TabularCapIntVar), sizeof(double)*num_var);
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
    td = scinew TypeDescription(TypeDescription::Other,
                                "TabularCapIntVar", true,
                                &makeMPI_TabularCapIntVar);
  }
  return td;
}
} // end namespace Uintah


