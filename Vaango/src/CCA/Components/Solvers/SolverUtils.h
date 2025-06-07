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

#ifndef __CCA_COMPONENTS_SOLVERS_SOLVERUTILS_H__
#define __CCA_COMPONENTS_SOLVERS_SOLVERUTILS_H__

#include <Core/Exceptions/ConvergenceFailure.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <Core/Disclosure/TypeUtils.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Grid/Variables/Array3.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/Stencil7.h>
#include <Core/Math/MiscMath.h>
#include <Core/Util/DebugStream.h>

namespace Uintah {

namespace Solver {

//  To turn on normal output
//  setenv SCI_DEBUG "CGSOLVER_UTIL_COUT:+"
inline Uintah::DebugStream cout_doing("CGSOLVER_UTIL_COUT", false);

void
Mult(Array3<double>& B,
     const Array3<Stencil7>& A,
     const Array3<double>& X,
     CellIterator iter,
     const IntVector& l,
     const IntVector& h1,
     long64& flops,
     long64& memrefs)
{
  if (cout_doing.active()) {
    cout_doing << "CGSolver::Mult" << std::endl;
  }

  for (; !iter.done(); ++iter) {
    IntVector idx      = *iter;
    const Stencil7& AA = A[idx];
    double result      = AA.p * X[idx];
    if (idx.x() > l.x()) {
      result += AA.w * X[idx + IntVector(-1, 0, 0)];
    }
    if (idx.x() < h1.x()) {
      result += AA.e * X[idx + IntVector(1, 0, 0)];
    }
    if (idx.y() > l.y()) {
      result += AA.s * X[idx + IntVector(0, -1, 0)];
    }
    if (idx.y() < h1.y()) {
      result += AA.n * X[idx + IntVector(0, 1, 0)];
    }
    if (idx.z() > l.z()) {
      result += AA.b * X[idx + IntVector(0, 0, -1)];
    }
    if (idx.z() < h1.z()) {
      result += AA.t * X[idx + IntVector(0, 0, 1)];
    }
    B[idx] = result;
  }
  IntVector diff = iter.end() - iter.begin();
  flops += 13 * diff.x() * diff.y() * diff.z();
  memrefs += 15L * diff.x() * diff.y() * diff.z() * 8L;
}

void
Mult(Array3<double>& B,
     const Array3<Stencil7>& A,
     const Array3<double>& X,
     CellIterator iter,
     const IntVector& l,
     const IntVector& h1,
     long64& flops,
     long64& memrefs,
     double& dotresult)
{
  if (cout_doing.active()) {
    cout_doing << "CGSolver::Mult" << std::endl;
  }

  double dot = 0;
  for (; !iter.done(); ++iter) {
    IntVector idx      = *iter;
    const Stencil7& AA = A[idx];
    double result      = AA.p * X[idx];
    if (idx.x() > l.x()) {
      result += AA.w * X[idx + IntVector(-1, 0, 0)];
    }
    if (idx.x() < h1.x()) {
      result += AA.e * X[idx + IntVector(1, 0, 0)];
    }
    if (idx.y() > l.y()) {
      result += AA.s * X[idx + IntVector(0, -1, 0)];
    }
    if (idx.y() < h1.y()) {
      result += AA.n * X[idx + IntVector(0, 1, 0)];
    }
    if (idx.z() > l.z()) {
      result += AA.b * X[idx + IntVector(0, 0, -1)];
    }
    if (idx.z() < h1.z()) {
      result += AA.t * X[idx + IntVector(0, 0, 1)];
    }
    B[idx] = result;
    dot += result * X[idx];
  }

  dotresult      = dot;
  IntVector diff = iter.end() - iter.begin();
  flops += 15 * diff.x() * diff.y() * diff.z();
  memrefs += 16L * diff.x() * diff.y() * diff.z() * 8L;
}

[[maybe_unused]] static void
Sub(Array3<double>& r,
    const Array3<double>& a,
    const Array3<double>& b,
    CellIterator iter,
    long64& flops,
    long64& memrefs)
{
  if (cout_doing.active()) {
    cout_doing << "CGSolver::Sub" << std::endl;
  }
  for (; !iter.done(); ++iter) {
    r[*iter] = a[*iter] - b[*iter];
  }
  IntVector diff = iter.end() - iter.begin();
  flops += diff.x() * diff.y() * diff.z();
  memrefs += diff.x() * diff.y() * diff.z() * 3L * 8L;
}

void
Mult(Array3<double>& r,
     const Array3<double>& a,
     const Array3<double>& b,
     CellIterator iter,
     long64& flops,
     long64& memrefs)
{
  if (cout_doing.active()) {
    cout_doing << "CGSolver::Mult" << std::endl;
  }
  for (; !iter.done(); ++iter) {
    r[*iter] = a[*iter] * b[*iter];
  }
  IntVector diff = iter.end() - iter.begin();
  flops += diff.x() * diff.y() * diff.z();
  memrefs += diff.x() * diff.y() * diff.z() * 3L * 8L;
}

[[maybe_unused]] static void
InverseDiagonal(Array3<double>& r,
                const Array3<Stencil7>& A,
                CellIterator iter,
                long64& flops,
                long64& memrefs)
{
  if (cout_doing.active()) {
    cout_doing << "CGSolver::InverseDiagonal" << std::endl;
  }
  for (; !iter.done(); ++iter) {
    r[*iter] = 1. / A[*iter].p;
  }
  IntVector diff = iter.end() - iter.begin();
  flops += diff.x() * diff.y() * diff.z();
  memrefs += 2L * diff.x() * diff.y() * diff.z() * 8L;
}

[[maybe_unused]] static double
L1(const Array3<double>& a, CellIterator iter, long64& flops, long64& memrefs)
{
  if (cout_doing.active()) {
    cout_doing << "CGSolver::L1" << std::endl;
  }
  double sum = 0;
  for (; !iter.done(); ++iter) {
    sum += Uintah::Abs(a[*iter]);
  }
  IntVector diff = iter.end() - iter.begin();
  flops += 2 * diff.x() * diff.y() * diff.z();
  memrefs += diff.x() * diff.y() * diff.z() * 8L;
  return sum;
}

double
LInf(const Array3<double>& a, CellIterator iter, long64& flops, long64& memrefs)
{
  if (cout_doing.active()) {
    cout_doing << "CGSolver::Linf" << std::endl;
  }
  double max = 0;
  for (; !iter.done(); ++iter) {
    max = Uintah::Max(max, Uintah::Abs(a[*iter]));
  }
  IntVector diff = iter.end() - iter.begin();
  flops += 2 * diff.x() * diff.y() * diff.z();
  memrefs += diff.x() * diff.y() * diff.z() * 8L;
  return max;
}

double
Dot(const Array3<double>& a,
    const Array3<double>& b,
    CellIterator iter,
    long64& flops,
    long64& memrefs)
{
  if (cout_doing.active()) {
    cout_doing << "CGSolver::Dot" << std::endl;
  }
  double sum = 0;
  for (; !iter.done(); ++iter) {
    sum += a[*iter] * b[*iter];
  }
  IntVector diff = iter.end() - iter.begin();
  flops += 2 * diff.x() * diff.y() * diff.z();
  memrefs += diff.x() * diff.y() * diff.z() * 2L * 8L;
  return sum;
}

void
ScMult_Add(Array3<double>& r,
           double s,
           const Array3<double>& a,
           const Array3<double>& b,
           CellIterator iter,
           long64& flops,
           long64& memrefs)
{
  if (cout_doing.active()) {
    cout_doing << "CGSolver::ScMult_Add" << std::endl;
  }

  for (; !iter.done(); ++iter) {
    r[*iter] = s * a[*iter] + b[*iter];
  }
  IntVector diff = iter.end() - iter.begin();
  flops += 2 * diff.x() * diff.y() * diff.z();
  memrefs += diff.x() * diff.y() * diff.z() * 3L * 8L;
}

} // namespace Solver
} // namespace Uintah

#endif //__CCA_COMPONENTS_SOLVERS_SOLVERUTILS_H__