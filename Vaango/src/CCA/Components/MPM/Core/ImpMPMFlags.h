/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __IMP_MPM_FLAGS_H__
#define __IMP_MPM_FLAGS_H__

#include <CCA/Components/MPM/Core/MPMFlags.h>

namespace Uintah {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class ImpMPMFlags
  \brief A structure that store the flags used for a MPM simulation
  \author John Schmidt \n
  C-SAFE and Department of Mechanical Engineering \n
  University of Utah \n
*/
/////////////////////////////////////////////////////////////////////////////

class ImpMPMFlags : public MPMFlags
{

public:
  ImpMPMFlags(const ProcessorGroup* myworld);

  virtual ~ImpMPMFlags() = default;

  virtual void readMPMFlags(ProblemSpecP& ps, Output* dataArchive);
  virtual void outputProblemSpec(ProblemSpecP& ps);

  bool d_projectHeatSource;
  bool d_doMechanics;
  double d_convCritDisp;
  double d_convCritEnergy;
  bool d_dynamic;
  int d_maxNumIterations;
  int d_numItersToDecreaseDelT;
  int d_numItersToIncreaseDelT;
  double d_delTDecreaseFactor;
  double d_delTIncreaseFactor;
  std::string d_solverType;
  bool d_tempSolve;
  bool d_interpolateParticleTempToGridEveryStep;

private:
  ImpMPMFlags(const ImpMPMFlags& state);
  ImpMPMFlags& operator=(const ImpMPMFlags& state);
};

} // End namespace Uintah

#endif // __IMP_MPM_FLAGS_H__
