/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 1997-2012 The University of Utah
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

/*
 *  MatrixOperations.h: Matrix Operations
 *
 *  Written by:
 *   Michael Callahan
 *   Department of Computer Science
 *   University of Utah
 *   August 2003
 *
 */

#ifndef __VAANGO_CORE_DATATYPE_MATRIX_OPERATIONS_H__
#define __VAANGO_CORE_DATATYPE_MATRIX_OPERATIONS_H__ 1

#include <Core/Datatypes/Matrix.h>

namespace Uintah {

MatrixHandle
operator+(MatrixHandle A, MatrixHandle B);
MatrixHandle
operator-(MatrixHandle A, MatrixHandle B);
MatrixHandle
operator*(MatrixHandle A, MatrixHandle B);
MatrixHandle
operator*(double a, MatrixHandle B);
MatrixHandle
operator*(MatrixHandle A, double b);

} // End namespace Uintah

#endif //__VAANGO_CORE_DATATYPE_MATRIX_OPERATIONS_H__
