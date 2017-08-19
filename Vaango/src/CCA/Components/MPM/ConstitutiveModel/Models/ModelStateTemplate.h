/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#ifndef __MODEL_STATE_TEMPLATE_H__
#define __MODEL_STATE_TEMPLATE_H__

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class ModelStateTemplate
  \brief Template class for a structure that stores the local state data.
  \author Biswajit Banerjee \n
*/
/////////////////////////////////////////////////////////////////////////////

template <class Model>
class ModelStateTemplate
{

public:
  Model d_model;

  ModelStateTemplate<Model>()
    : d_model()
  {
  }

  ModelStateTemplate<Model>(const ModelStateTemplate<Model>& state)
  {
    d_model = state.d_model;
  }

  ModelStateTemplate<Model>(const ModelStateTemplate<Model>* state)
  {
    d_model = state->d_model;
  }

  ~ModelStateTemplate<Model>() {}

  ModelStateTemplate<Model>& operator=(const ModelStateTemplate<Model>& state)
  {
    if (this == &state)
      return *this;
    d_model = state.d_model;
    return *this;
  }

  ModelStateTemplate<Model>* operator=(const ModelStateTemplate<Model>* state)
  {
    if (this == state)
      return *this;
    d_model = state->d_model;
    return this;
  }
};

} // End namespace Vaango

#endif // __MODEL_STATE_TEMPLATE_H__
