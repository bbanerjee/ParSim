/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <Core/Disclosure/TypeUtils.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Variables/ReductionVariable.h>
#include <Core/Grid/Variables/Reductions.h>
#include <Core/Util/FancyAssert.h>

#include <sci_defs/bits_defs.h> // for SCI_32BITS
#include <sci_defs/osx_defs.h>  // for OSX_SNOW_LEOPARD_OR_LATER

namespace Uintah {

template<>
void
ReductionVariable<double, Reductions::Min<double>>::getMPIInfo(
  int& count,
  MPI_Datatype& datatype,
  MPI_Op& op)
{
  datatype = MPI_DOUBLE;
  count    = 1;
  op       = MPI_MIN;
}

template<>
void
ReductionVariable<double, Reductions::Min<double>>::getMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(double)));
  double* ptr = reinterpret_cast<double*>(&data[index]);
  *ptr        = *(d_value.get());
  index += sizeof(double);
}

#if !defined(SCI_32BITS)
template<>
void
ReductionVariable<long long, Reductions::Min<long long>>::getMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(long long)));
  long long* ptr = reinterpret_cast<long long*>(&data[index]);
  *ptr           = *(d_value.get());
  index += sizeof(long long);
}

template<>
void
ReductionVariable<long long, Reductions::Sum<long long>>::getMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(long long)));
  long long* ptr = reinterpret_cast<long long*>(&data[index]);
  *ptr           = *(d_value.get());
  index += sizeof(long long);
}
#endif

template<>
void
ReductionVariable<double, Reductions::Min<double>>::putMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(double)));
  double* ptr = reinterpret_cast<double*>(&data[index]);
  d_value     = std::make_shared<double>(*ptr);
  index += sizeof(double);
}

template<>
void
ReductionVariable<double, Reductions::Max<double>>::getMPIInfo(
  int& count,
  MPI_Datatype& datatype,
  MPI_Op& op)
{
  datatype = MPI_DOUBLE;
  count    = 1;
  op       = MPI_MAX;
}

template<>
void
ReductionVariable<double, Reductions::Max<double>>::getMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(double)));
  double* ptr = reinterpret_cast<double*>(&data[index]);
  *ptr        = *(d_value.get());
  index += sizeof(double);
}

template<>
void
ReductionVariable<double, Reductions::Max<double>>::putMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(double)));
  double* ptr = reinterpret_cast<double*>(&data[index]);
  d_value     = std::make_shared<double>(*ptr);
  index += sizeof(double);
}

template<>
void
ReductionVariable<double, Reductions::Sum<double>>::getMPIInfo(
  int& count,
  MPI_Datatype& datatype,
  MPI_Op& op)
{
  datatype = MPI_DOUBLE;
  count    = 1;
  op       = MPI_SUM;
}

#if !defined(SCI_32BITS)
template<>
void
ReductionVariable<long long, Reductions::Sum<long long>>::getMPIInfo(
  int& count,
  MPI_Datatype& datatype,
  MPI_Op& op)
{
  datatype = MPI_LONG_LONG;
  count    = 1;
  op       = MPI_SUM;
}
#endif

template<>
void
ReductionVariable<double, Reductions::Sum<double>>::getMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(double)));
  double* ptr = reinterpret_cast<double*>(&data[index]);
  *ptr        = *(d_value.get());
  index += sizeof(double);
}

template<>
void
ReductionVariable<double, Reductions::Sum<double>>::putMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(double)));
  double* ptr = reinterpret_cast<double*>(&data[index]);
  d_value     = std::make_shared<double>(*ptr);
  index += sizeof(double);
}

template<>
void
ReductionVariable<bool, Reductions::And<bool>>::getMPIInfo(
  int& count,
  MPI_Datatype& datatype,
  MPI_Op& op)
{
  datatype = MPI_CHAR;
  count    = 1;
  op       = MPI_LAND;
}

template<>
void
ReductionVariable<bool, Reductions::And<bool>>::getMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(char)));
  char* ptr = reinterpret_cast<char*>(&data[index]);
  *ptr      = *(d_value.get());
  index += sizeof(char);
}

template<>
void
ReductionVariable<bool, Reductions::And<bool>>::putMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(char)));
  char* ptr = reinterpret_cast<char*>(&data[index]);
  d_value   = std::make_shared<bool>(*ptr);
  index += sizeof(char);
}

template<>
void
ReductionVariable<bool, Reductions::Or<bool>>::getMPIInfo(
  int& count,
  MPI_Datatype& datatype,
  MPI_Op& op)
{
  datatype = MPI_CHAR;
  count    = 1;
  op       = MPI_LOR;
}

template<>
void
ReductionVariable<bool, Reductions::Or<bool>>::getMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(char)));
  char* ptr = reinterpret_cast<char*>(&data[index]);
  *ptr      = *(d_value.get());
  index += sizeof(char);
}

template<>
void
ReductionVariable<bool, Reductions::Or<bool>>::putMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(char)));
  char* ptr = reinterpret_cast<char*>(&data[index]);
  d_value     = std::make_shared<bool>(*ptr);
  index += sizeof(char);
}

template<>
void
ReductionVariable<long64, Reductions::Sum<long64>>::getMPIInfo(
  int& count,
  MPI_Datatype& datatype,
  MPI_Op& op)
{
  datatype = MPI_LONG;
  count    = 1;
  op       = MPI_SUM;
}

template<>
void
ReductionVariable<long64, Reductions::Sum<long64>>::getMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(long)));
  long* ptr = reinterpret_cast<long*>(&data[index]);
  *ptr      = *(d_value.get());
  index += sizeof(long);
}

template<>
void
ReductionVariable<long64, Reductions::Sum<long64>>::putMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(long)));
  long* ptr = reinterpret_cast<long*>(&data[index]);
  d_value   = std::make_shared<long64>(*ptr);
  index += sizeof(long);
}

template<>
void
ReductionVariable<long long, Reductions::Sum<long long>>::putMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - sizeof(long long)));
  long long* ptr = reinterpret_cast<long long*>(&data[index]);
  d_value        = std::make_shared<long long>(*ptr);
  index += sizeof(long long);
}

template<>
void
ReductionVariable<Vector, Reductions::Sum<Vector>>::getMPIInfo(
  int& count,
  MPI_Datatype& datatype,
  MPI_Op& op)
{
  datatype = MPI_DOUBLE;
  count    = 3;
  op       = MPI_SUM;
}

template<>
void
ReductionVariable<Vector, Reductions::Sum<Vector>>::getMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - 3 * sizeof(double)));
  double* ptr = reinterpret_cast<double*>(&data[index]);
  *ptr++      = d_value.get()->x();
  *ptr++      = d_value.get()->y();
  *ptr++      = d_value.get()->z();
  index += 3 * sizeof(double);
}

template<>
void
ReductionVariable<Vector, Reductions::Sum<Vector>>::putMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - 3 * sizeof(double)));
  double* ptr = reinterpret_cast<double*>(&data[index]);
  d_value.get()->x(*ptr++);
  d_value.get()->y(*ptr++);
  d_value.get()->z(*ptr++);
  index += 3 * sizeof(double);
}

template<>
void
ReductionVariable<Vector, Reductions::Min<Vector>>::getMPIInfo(
  int& count,
  MPI_Datatype& datatype,
  MPI_Op& op)
{
  datatype = MPI_DOUBLE;
  count    = 3;
  op       = MPI_MIN;
}

template<>
void
ReductionVariable<Vector, Reductions::Min<Vector>>::getMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - 3 * sizeof(double)));
  double* ptr = reinterpret_cast<double*>(&data[index]);
  *ptr++      = d_value.get()->x();
  *ptr++      = d_value.get()->y();
  *ptr++      = d_value.get()->z();
  index += 3 * sizeof(double);
}

template<>
void
ReductionVariable<Vector, Reductions::Min<Vector>>::putMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - 3 * sizeof(double)));
  double* ptr = reinterpret_cast<double*>(&data[index]);
  d_value.get()->x(*ptr++);
  d_value.get()->y(*ptr++);
  d_value.get()->z(*ptr++);
  index += 3 * sizeof(double);
}

template<>
void
ReductionVariable<Vector, Reductions::Max<Vector>>::getMPIInfo(
  int& count,
  MPI_Datatype& datatype,
  MPI_Op& op)
{
  datatype = MPI_DOUBLE;
  count    = 3;
  op       = MPI_MAX;
}

template<>
void
ReductionVariable<Vector, Reductions::Max<Vector>>::getMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - 3 * sizeof(double)));
  double* ptr = reinterpret_cast<double*>(&data[index]);
  *ptr++      = d_value.get()->x();
  *ptr++      = d_value.get()->y();
  *ptr++      = d_value.get()->z();
  index += 3 * sizeof(double);
}

template<>
void
ReductionVariable<Vector, Reductions::Max<Vector>>::putMPIData(
  std::vector<char>& data,
  int& index)
{
  ASSERTRANGE(index, 0, static_cast<int>(data.size() + 1 - 3 * sizeof(double)));
  double* ptr = reinterpret_cast<double*>(&data[index]);
  d_value.get()->x(*ptr++);
  d_value.get()->y(*ptr++);
  d_value.get()->z(*ptr++);
  index += 3 * sizeof(double);
}

} // end namespace Uintah
