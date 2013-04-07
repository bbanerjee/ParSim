#ifndef EMU2DC_VECTOR3_H
#define EMU2DC_VECTOR3_H

#include <Types.h>

namespace Emu2DC {
  class Vector3 : public Array3 {
    public:
    
    Vector3 operator+(const Vector3 vec) {
      Vector3 out_vec;
      out_vec[0] = (*this)[0] + vec[0]; 
      out_vec[1] = (*this)[1] + vec[1]; 
      out_vec[2] = (*this)[2] + vec[2]; 
      return out_vec;
    }
  };
}

#endif

