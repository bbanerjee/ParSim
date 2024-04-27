/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#ifndef __DEFORMATION_STATE_DATA_H__
#define __DEFORMATION_STATE_DATA_H__

#include <Core/Math/Matrix3.h>

namespace Uintah {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class DeformationState
  \brief A structure that stores deformation state data
  \author Todd Harman \n
*/
/////////////////////////////////////////////////////////////////////////////

struct DeformationState
{
  Matrix3 D{ 0.0 };    // Rate of deformation
  Matrix3 devD{ 0.0 }; // Deviatoric part of rate of deformation
  Matrix3 devStressInc{ 0.0 };
  double J{ 1.0 };
  double viscoElasticWorkRate{ 0.0 };
  ~DeformationState() = default;
};
} // End namespace Uintah

#endif // __DEFORMATION_STATE_DATA_H__
