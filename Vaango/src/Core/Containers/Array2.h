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

/*
 *  Array2.h: Interface to dynamic 2D array class
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 */

#ifndef __CORE_CONTAINERS_ARRAY2_H__
#define __CORE_CONTAINERS_ARRAY2_H__ 1

#include <Core/Util/Assert.h>

namespace Uintah {

class RigorousTest;

template<class T>
class Array2
{
private:
  T** objs;
  int dm1;
  int dm2;
  void
  allocate();

  // The copy constructor and the assignment operator have been
  // privatized on purpose -- no one should use these.  Instead,
  // use the default constructor and the copy method.
  //////////
  // Array2 Copy Constructor
  Array2(const Array2&);
  //////////
  // Assignment Operator
  Array2<T>&
  operator=(const Array2&);

public:
  //////////
  // Create a 0X0 Array
  Array2();

  //////////
  // Create an n by n array
  Array2(int, int);

  //////////
  // Class Destructor
  virtual ~Array2();

  //////////
  // Used for accessing elements in the Array
  inline T&
  operator()(int d1, int d2) const
  {
    ASSERTL3(d1 >= 0 && d1 < dm1);
    ASSERTL3(d2 >= 0 && d2 < dm2);
    return objs[d1][d2];
  }

  //////////
  // Array2 Copy Method
  void
  copy(const Array2&);

  //////////
  // Returns number of rows
  inline int
  dim1() const
  {
    return dm1;
  }

  //////////
  // Returns number of cols
  inline int
  dim2() const
  {
    return dm2;
  }

  //////////
  // Resize Array
  void
  resize(int, int);

  //////////
  // Initialize all values in an array
  void
  initialize(const T&);

  inline T**
  get_dataptr()
  {
    return objs;
  }

};

template<class T>
Array2<T>::Array2()
{
  dm1 = dm2 = 0;
  objs      = 0;
}

template<class T>
void
Array2<T>::allocate()
{
  if (dm1 == 0 || dm2 == 0) {
    objs = 0;
  } else {
    objs = new T*[dm1];
    T* p = new T[dm1 * dm2];
    for (int i = 0; i < dm1; i++) {
      objs[i] = p;
      p += dm2;
    }
  }
}

template<class T>
void
Array2<T>::resize(int d1, int d2)
{
  if (objs && dm1 == d1 && dm2 == d2) {
    return;
  }
  dm1 = d1;
  dm2 = d2;
  if (objs) {
    delete[] objs[0];
    delete[] objs;
  }
  allocate();
}

template<class T>
Array2<T>::Array2(const Array2<T>& a)
  : dm1(a.dm1)
  , dm2(a.dm2)
{
  allocate();
}

template<class T>
Array2<T>::Array2(int dm1, int dm2)
  : dm1(dm1)
  , dm2(dm2)
{
  allocate();
}

template<class T>
Array2<T>::~Array2()
{
  if (objs) {
    delete[] objs[0];
    delete[] objs;
  }
}

template<class T>
void
Array2<T>::initialize(const T& t)
{
  ASSERT(dm1 == 0 || dm2 == 0 || objs != 0);
  for (int i = 0; i < dm1; i++) {
    for (int j = 0; j < dm2; j++) {
      objs[i][j] = t;
    }
  }
}

template<class T>
void
Array2<T>::copy(const Array2<T>& copy)
{
  resize(copy.dim1(), copy.dim2());
  for (int i = 0; i < dm1; i++) {
    for (int j = 0; j < dm2; j++) {
      objs[i][j] = copy.objs[i][j];
    }
  }
}


} // End namespace Uintah

#endif //__CORE_CONTAINERS_ARRAY2_H__
