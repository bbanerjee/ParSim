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

#include <ShapeFunctions/MPMShapeFunctionFactory.h>
#include <Exception.h>

#include <ShapeFunctions/MPMShapeFunction.h>
#include <ShapeFunctions/LinearShapeFunction.h>
#include <ShapeFunctions/GIMPShapeFunction.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <string>

using namespace BrMPM;

MPMShapeFunctionFactory::MPMShapeFunctionFactory()
{
}

MPMShapeFunctionFactory::~MPMShapeFunctionFactory()
{
}

MPMShapeFunctionP
MPMShapeFunctionFactory::create(const Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP shape_ps = ps->findBlock("ShapeFunction");
  if (!shape_ps) {
    throw Exception("**ERROR** <ShapeFunction> tag not found", __FILE__, __LINE__);
  }
 
  std::string shape_type;
  if (!shape_ps->getAttribute("type", shape_type)) {
    throw Exception("**ERROR** <ShapeFunction type=?> Type attribute not found",
                    __FILE__, __LINE__);
  }

  // Create shape function
  if (shape_type == "GIMP") {
    return std::make_shared<GIMPShapeFunction>();
  } else if (shape_type == "Linear") {
    return std::make_shared<LinearShapeFunction>();
  } else {
    std::ostringstream out;
    out << "**ERROR** Unknown shape function" << shape_type << " for MPM. " << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }
}


      















  
