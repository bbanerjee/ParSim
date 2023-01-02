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


#include <Core/Geometry/IntVector.h>
#include <Core/Persistent/Persistent.h>
#include <Core/Util/FETypeDescription.h>
#include <Core/Util/XMLUtils.h>

#include <stdlib.h>
#include <iostream>
using std::ostream;

namespace Uintah{

  void
  Pio(Piostream& stream, IntVector& p)
  {
    stream.begin_cheap_delim();
    Pio(stream, p.value_[0]);
    Pio(stream, p.value_[1]);
    Pio(stream, p.value_[2]);
    stream.end_cheap_delim();
  }



  const string& 
  IntVector::get_h_file_path() {
    static const string path(FETypeDescription::cc_to_h(__FILE__));
    return path;
  }

  const FETypeDescription*
  get_fetype_description(IntVector*)
  {
    static FETypeDescription* td = 0;
    if(!td){
      td = scinew FETypeDescription("IntVector", IntVector::get_h_file_path(), "Uintah");
    }
    return td;
  }

  ostream&
  operator<<(std::ostream& out, const Uintah::IntVector& v)
  {
    out << "[int " << v.x() << ", " << v.y() << ", " << v.z() << ']';
    return out;
  }

  IntVector
  IntVector::fromString( const string & source )
  {
    IntVector result;

    // Parse out the [num,num,num]

    string::size_type i1 = source.find("[");
    string::size_type i2 = source.find_first_of(",");
    string::size_type i3 = source.find_last_of(",");
    string::size_type i4 = source.find("]");
  
    string x_val(source,i1+1,i2-i1-1);
    string y_val(source,i2+1,i3-i2-1);
    string z_val(source,i3+1,i4-i3-1);

    validateType( x_val, Uintah::UintahXML::INT_TYPE );
    validateType( y_val, Uintah::UintahXML::INT_TYPE );
    validateType( z_val, Uintah::UintahXML::INT_TYPE );
          
    result.x( atoi(x_val.c_str()) );
    result.y( atoi(y_val.c_str()) );
    result.z( atoi(z_val.c_str()) );

    return result;
  }

} //end namespace Uintah

