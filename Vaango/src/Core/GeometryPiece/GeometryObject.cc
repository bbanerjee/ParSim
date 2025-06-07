/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

namespace Uintah {

GeometryObject::GeometryObject(GeometryPieceP piece,
                               ProblemSpecP& ps,
                               std::list<DataItem>& data)
    : d_piece(piece) {
  for (const auto& dataItem : data) {
    // std::cout << "\t Read " << dataItem.type << " name = " << dataItem.name
    // << std::endl;
    switch (dataItem.type) {
      case Double: {
        double val;
        if (dataItem.name == "volumeFraction") {
          ps->getWithDefault(dataItem.name, val, -1.0);
        } else {
          ps->require(dataItem.name, val);
        }
        d_double_data[dataItem.name] = val;
        break;
      }
      case Integer: {
        int val;
        ps->require(dataItem.name, val);
        d_int_data[dataItem.name] = val;
        break;
      }
      case Vector: {
        Uintah::Vector val;
        if (dataItem.name == "affineTransformation_A0") {
          ps->getWithDefault(dataItem.name, val, Uintah::Vector(1., 0., 0.));

        } else if (dataItem.name == "affineTransformation_A1") {
          ps->getWithDefault(dataItem.name, val, Uintah::Vector(0., 1., 0.));

        } else if (dataItem.name == "affineTransformation_A2") {
          ps->getWithDefault(dataItem.name, val, Uintah::Vector(0., 0., 1.));

        } else if (dataItem.name == "affineTransformation_b") {
          ps->getWithDefault(dataItem.name, val, Uintah::Vector(0., 0., 0.));

        } else {
          ps->require(dataItem.name, val);
        }
        d_vector_data[dataItem.name] = val;
        break;
      }
      case IntVector: {
        Uintah::IntVector val;
        ps->require(dataItem.name, val);
        d_intvector_data[dataItem.name] = val;
        // std::cout << "\t Read IntVector" << dataItem.name << " = " << val <<
        // std::endl;
        break;
      }
      case Point: {
        Uintah::Point val;
        ps->require(dataItem.name, val);
        d_point_data[dataItem.name] = val;
        break;
      }
    };
  }
}

void
GeometryObject::outputProblemSpec(ProblemSpecP& ps) {
  ProblemSpecP geom_obj_ps = ps->appendChild("geom_object");
  d_piece->outputProblemSpec(geom_obj_ps);

  for (const auto& data : d_double_data) {
    if (!(data.first.compare("volumeFraction") == 0 && data.second == -1.0))
      geom_obj_ps->appendElement(data.first.c_str(), data.second);
  }
  for (const auto& data : d_vector_data) {
    geom_obj_ps->appendElement(data.first.c_str(), data.second);
  }
  for (const auto& data : d_intvector_data) {
    geom_obj_ps->appendElement(data.first.c_str(), data.second);
  }
  for (const auto& data : d_point_data) {
    geom_obj_ps->appendElement(data.first.c_str(), data.second);
  }
}

}  // end namespace Uintah