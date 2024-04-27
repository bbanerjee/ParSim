/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 1997-2012 The University of Utah
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

#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/BoundaryConditions/BoundCond.h>
#include <Core/Grid/BoundaryConditions/BoundCondFactory.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

namespace Uintah {

auto
BoundCondFactory::create(ProblemSpecP& child,
                         int& mat_id,
                         const std::string face_label) -> BoundCondBaseSP

{
  std::map<std::string, std::string> bc_attr;
  child->getAttributes(bc_attr);

  // Check to see if "id" is defined
  if (bc_attr.find("id") == bc_attr.end()) {
    SCI_THROW(ProblemSetupException("id is not specified in the BCType tag",
                                    __FILE__,
                                    __LINE__));
  }

  if (bc_attr["id"] != "all") {
    std::istringstream ss(bc_attr["id"]);
    ss >> mat_id;
  } else {
    mat_id = -1;
  }

  //  std::cout << "mat_id = " << mat_id << std::endl;
  // Determine whether or not things are a scalar, Vector or a NoValue, i.e.
  // Symmetry

  int i_value;
  double d_value;
  Vector v_value;
  std::string s_value = "none";

  std::string valAttribute;
  bool attrPS          = child->getAttribute("value", valAttribute);
  ProblemSpecP valuePS = child->findBlock("value");

  if (valuePS && attrPS) {
    // the user specified BOTH <value> </value> AND value="" attribute
    SCI_THROW(ProblemSetupException(
      "Error: It looks like you specified two values for BC " +
        bc_attr["label"] +
        ". This is not allowed! You could only specify either a value "
        "attribute or a <value> node. Please revise your input file.",
      __FILE__,
      __LINE__));
  }

  BoundCondBaseSP bc;

  if (attrPS) {
    ProblemSpec::InputType theInputType = child->getInputType(valAttribute);
    switch (theInputType) {
      case ProblemSpec::NUMBER_TYPE: {
        if (bc_attr["type"] == "int") {
          // integer ONLY if the tag 'type = "int"' is added
          child->getAttribute("value", i_value);
          bc = std::make_shared<BoundCond<int>>(bc_attr["label"],
                                                bc_attr["var"],
                                                i_value,
                                                face_label,
                                                BoundCondBase::INT_TYPE);
        } else { // double (default)
          child->getAttribute("value", d_value);
          bc = std::make_shared<BoundCond<double>>(bc_attr["label"],
                                                   bc_attr["var"],
                                                   d_value,
                                                   face_label,
                                                   BoundCondBase::DOUBLE_TYPE);
        }
      } break;
      case ProblemSpec::VECTOR_TYPE:
        child->getAttribute("value", v_value);
        bc = std::make_shared<BoundCond<Vector>>(bc_attr["label"],
                                                 bc_attr["var"],
                                                 v_value,
                                                 face_label,
                                                 BoundCondBase::VECTOR_TYPE);
        break;
      case ProblemSpec::STRING_TYPE:
        bc =
          std::make_shared<BoundCond<std::string>>(bc_attr["label"],
                                                   bc_attr["var"],
                                                   valAttribute,
                                                   face_label,
                                                   BoundCondBase::STRING_TYPE);
        break;
      case ProblemSpec::UNKNOWN_TYPE:
      default:
        bc = std::make_shared<BoundCond<NoValue>>(bc_attr["label"],
                                                  bc_attr["var"]);
        break;
    }
  } else if (valuePS) { // Found <value> tag.
    child->get("value", s_value);
    ProblemSpec::InputType theInputType = child->getInputType(s_value);

    switch (theInputType) {
      case ProblemSpec::NUMBER_TYPE:
        if (bc_attr["type"] == "int") {
          child->get("value", i_value);
          bc = std::make_shared<BoundCond<double>>(bc_attr["label"],
                                                   bc_attr["var"],
                                                   i_value,
                                                   face_label,
                                                   BoundCondBase::INT_TYPE);
        } else {
          child->get("value", d_value);
          bc = std::make_shared<BoundCond<double>>(bc_attr["label"],
                                                   bc_attr["var"],
                                                   d_value,
                                                   face_label,
                                                   BoundCondBase::DOUBLE_TYPE);
        }
        break;
      case ProblemSpec::VECTOR_TYPE:
        child->get("value", v_value);
        bc = std::make_shared<BoundCond<Vector>>(bc_attr["label"],
                                                 bc_attr["var"],
                                                 v_value,
                                                 face_label,
                                                 BoundCondBase::VECTOR_TYPE);
        break;
      case ProblemSpec::STRING_TYPE:
        bc =
          std::make_shared<BoundCond<std::string>>(bc_attr["label"],
                                                   bc_attr["var"],
                                                   s_value,
                                                   face_label,
                                                   BoundCondBase::STRING_TYPE);
        break;
      case ProblemSpec::UNKNOWN_TYPE:
      default:
        bc = std::make_shared<BoundCond<NoValue>>(bc_attr["label"],
                                                  bc_attr["var"]);
        break;
    }
  } else {
    bc = std::make_shared<BoundCond<NoValue>>(bc_attr["label"], bc_attr["var"]);
  }
  if (!bc) {
    throw InternalError("Problem creating boundary condition",
                        __FILE__,
                        __LINE__);
  }
  return bc;
}

auto
BoundCondFactory::customBC([[maybe_unused]] int mat_id,
                           const std::string face_label,
                           double value,
                           const std::string label,
                           const std::string var) -> BoundCondBaseSP
{
  return std::make_shared<BoundCond<double>>(label,
                                             var,
                                             value,
                                             face_label,
                                             BoundCondBase::DOUBLE_TYPE);
}

auto
BoundCondFactory::customBC([[maybe_unused]] int mat_id,
                           const std::string face_label,
                           const Vector value,
                           const std::string label,
                           const std::string var) -> BoundCondBaseSP
{
  return std::make_shared<BoundCond<Vector>>(label,
                                             var,
                                             value,
                                             face_label,
                                             BoundCondBase::VECTOR_TYPE);
}

auto
BoundCondFactory::customBC([[maybe_unused]] int mat_id,
                           const std::string face_label,
                           const std::string value,
                           const std::string label,
                           const std::string var) -> BoundCondBaseSP
{
  return std::make_shared<BoundCond<std::string>>(label,
                                                  var,
                                                  value,
                                                  face_label,
                                                  BoundCondBase::STRING_TYPE);
}

} // namespace Uintah