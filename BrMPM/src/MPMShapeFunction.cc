#include <MPMShapeFunction.h>
#include <Exception.h>
#include <iostream>
#include <string>

#include <Core/ProblemSpec/ProblemSpec.h>

using namespace MPM;

MPMShapeFunction::MPMShapeFunction()
            : d_shape_size(9), d_ghost(2)
{
}

MPMShapeFunction::~MPMShapeFunction()
{
}

void
MPMShapeFunction::initialise(const Uintah::ProblemSpecP& ps)
{
  Uintah::problemSpecP shape_ps = ps->findblock("ShapeFunction");
  if (!shape_ps) {
    throw Exception("**ERROR** <ShapeFunction> tag not found", __FILE__, __LINE__);
  }
 
  std::string shape_type;
  shape_ps->require("shape_function", shape_type);
  if (shape_type == "GIMP") {
    d_shape = ShapeType::GIMP;
    d_shape_size = 9;
    d_ghost = 2;
  } else if (shape_type == "Quad") {
    d_shape = ShapeType::Quad;
    d_shape_size = 9;
    d_ghost = 2;
  } else if (shape_type == "Linear") {
    d_shape = ShapeType::Linear;
    d_shape_size = 4;
    d_ghost = 1;
  } else if (shape_type == "Cubic") {
    d_shape = ShapeType::Cubic;
    d_shape_size = 12;
    d_ghost = 2;
  } else {
    std::ostringstream out;
    out << "**ERROR** UNknown shape function" << shape_type << " for MPM. " << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

}


      















  
