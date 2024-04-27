/*
 * The MIT License
 *
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef __MPM_INTERNAL_VAR_TYPES_H__
#define __MPM_INTERNAL_VAR_TYPES_H__

#include <Core/Disclosure/TypeUtils.h>
#include <Core/Math/Matrix3.h>
#include <Core/Util/Assert.h>

namespace Uintah {

class TypeDescription;

class MetalIntVar
{
public:
  double eqPlasticStrain;
  double plasticPorosity;
};

class DStressDMetalIntVar
{
public:
  Matrix3 eqPlasticStrain;
  Matrix3 plasticPorosity;
};

class ArenaIntVar
{
public:
  double kappa;
  double capX;
  double plasticVolStrain;
  double P3;
  Matrix3 plasticStrain;
};

class BorjaIntVar
{
public:
  double p_c;
};

class SoilBrannonIntVar
{
public:
  double kappa;
};

class TabularCapIntVar
{
public:
  double capX;
};

struct ViscoScramStateData
{
  Matrix3 DevStress[5];
};

} // End namespace Uintah

#include <Core/Datatypes/TypeName.h>
namespace Uintah {

void
swapbytes(Uintah::MetalIntVar& mp);

template<>
const std::string
find_type_name(Uintah::MetalIntVar* mp);

const TypeDescription*
fun_getTypeDescription(Uintah::MetalIntVar*);

void
swapbytes(Uintah::DStressDMetalIntVar& mp);

template<>
const std::string
find_type_name(Uintah::DStressDMetalIntVar* mp);

const TypeDescription*
fun_getTypeDescription(Uintah::DStressDMetalIntVar*);

void
swapbytes(Uintah::ArenaIntVar& mp);
template<>

const std::string
find_type_name(Uintah::ArenaIntVar* mp);

const TypeDescription*
fun_getTypeDescription(Uintah::ArenaIntVar*);

void
swapbytes(Uintah::BorjaIntVar& mp);
template<>

const std::string
find_type_name(Uintah::BorjaIntVar* mp);

const TypeDescription*
fun_getTypeDescription(Uintah::BorjaIntVar*);

void
swapbytes(Uintah::SoilBrannonIntVar& mp);

template<>
const std::string
find_type_name(Uintah::SoilBrannonIntVar* mp);

const TypeDescription*
fun_getTypeDescription(Uintah::SoilBrannonIntVar*);

void
swapbytes(Uintah::TabularCapIntVar& mp);

template<>
const std::string
find_type_name(Uintah::TabularCapIntVar* mp);

const TypeDescription*
fun_getTypeDescription(Uintah::TabularCapIntVar*);

/*! Set up type for StateData */
void
swapbytes(Uintah::ViscoScramStateData& d);

template<>
const std::string
find_type_name(Uintah::ViscoScramStateData* mp);

const TypeDescription*
fun_getTypeDescription(Uintah::ViscoScramStateData*);

} // End namespace Uintah

#endif //__MPM_INTERNAL_VAR_TYPES_H__
