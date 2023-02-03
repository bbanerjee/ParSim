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

#ifndef __VAANGO_PERIDYNAMICS_FLAGS_H__
#define __VAANGO_PERIDYNAMICS_FLAGS_H__

#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <CCA/Ports/Output.h>
#include <Core/Geometry/Vector.h>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class PeridynamicsFlags
    \brief A structure that store the flags used for a Peridynamics simulation
  */
  /////////////////////////////////////////////////////////////////////////////


  class PeridynamicsFlags {

  public:

    enum IntegratorType {
      ForwardEuler,
      VelocityVerlet,
      BackwardEuler,
      None
    };

    const Uintah::ProcessorGroup* d_myworld;

    Uintah::Vector d_gravity;
    std::string d_integratorType; // Explicit or implicit time integration
    IntegratorType d_integrator;
    double d_numCellsInHorizon;
    bool d_useLoadCurves;

    PeridynamicsFlags(const Uintah::ProcessorGroup* myworld);

    virtual ~PeridynamicsFlags();

    virtual void readPeridynamicsFlags(Uintah::ProblemSpecP& ps, Uintah::Output* dataArchive);
    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);

  private:

    PeridynamicsFlags(const PeridynamicsFlags& state);
    PeridynamicsFlags& operator=(const PeridynamicsFlags& state);
    
  };

} // End namespace Vaango

#endif  // __VAANGO_PERIDYNAMICS_FLAGS_H__ 

