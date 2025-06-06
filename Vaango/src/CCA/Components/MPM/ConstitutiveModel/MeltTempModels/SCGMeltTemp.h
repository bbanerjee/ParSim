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

/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#ifndef __SCG_MELT_TEMP_MODEL_H__
#define __SCG_MELT_TEMP_MODEL_H__

#include "MeltingTempModel.h"
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

/*! \class SCGMeltTemp
 *  \brief The melt temp model used by Steinberg,Cochran,Guinan in
 *         the SCG plasticity model.
 *  \author Biswajit Banerjee,
 *  \author C-SAFE and Department of Mechanical Engineering,
 *  \author University of Utah.
 *
*/
class SCGMeltTemp : public MeltingTempModel
{

private:
  double d_Gamma0; // Material constant (also in SCG model)
  double d_a;      // Material constant (also in SCG model)
  double d_Tm0;    // Material constant (also in SCG model)

  SCGMeltTemp& operator=(const SCGMeltTemp& mtm);

public:
  /*! Construct a constant melt temp model. */
  SCGMeltTemp(ProblemSpecP& ps);

  /*! Construct a copy of constant melt temp model. */
  SCGMeltTemp(const SCGMeltTemp* mtm);

  /*! Destructor of constant melt temp model.   */
  ~SCGMeltTemp() override;

  void outputProblemSpec(ProblemSpecP& ps) override;

  /*! Compute the melt temp */
  double computeMeltingTemp(const ModelStateBase* state) override;
};
} // End namespace Uintah

#endif // __SCG_MELT_TEMP_MODEL_H__
