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
 *  Array1.h: Interface to dynamic 1D array class
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 */

#ifndef __CORE_CONTAINERS_ARRAY1_H__
#define __CORE_CONTAINERS_ARRAY1_H__ 1

#include <Core/Util/Assert.h>

namespace Uintah {

class RigorousTest;

template<class T>
class Array1
{
  T* objs;
  int _size;
  int nalloc;
  int default_grow_size;

public:
  //////////
  // Copy the array - this can be costly, so try to avoid it.
  Array1(const Array1&);

  //////////
  // Make a new array 1. <i>size</i> gives the initial size of the array,
  //<i>default_grow_size</i> indicates the minimum number of objects that
  // should be added to the array at a time.  <i>asize</i> tells how many
  // objects should be allocated initially
  Array1(int size = 0, int default_grow_size = 10, int asize = -1);

  //////////
  // Copy over the array - this can be costly, so try to avoid it.
  Array1<T>&
  operator=(const Array1&);

  //////////
  // Compare over the array - this can be costly, so try to avoid it.
  int
  operator==(const Array1<T>&) const;
  int
  operator!=(const Array1<T>&) const;

  //////////
  // Deletes the array and frees the associated memory
  ~Array1();

  //////////
  // Accesses the nth element of the array
  inline const T&
  operator[](int n) const
  {
    CHECKARRAYBOUNDS(n, 0, _size);
    return objs[n];
  }

  //////////
  // Accesses the nth element of the array
  inline T&
  operator[](int n)
  {
    CHECKARRAYBOUNDS(n, 0, _size);
    return objs[n];
  }

  //////////
  // Returns the size of the array
  inline int
  size() const
  {
    return _size;
  }

  //////////
  // Make the array larger by count elements
  void
  grow(int count, int grow_size = 10);

  //////////
  // set allocated size
  void
  reserve(int n);

  //////////
  // Add one element to the array.  equivalent to:
  //  grow(1)
  //  array[array.size()-1]=data
  void
  add(const T&);

  //////////
  // Insert one element into the array.  This is very inefficient
  // if you insert anything besides the (new) last element.
  void
  insert(int, const T&);

  //////////
  // Remove one element from the array.  This is very inefficient
  // if you remove anything besides the last element.
  void
  remove(int);

  //////////
  // Remove all elements in the array.  The array is not freed,
  // and the number of allocated elements remains the same.
  void
  remove_all();

  //////////
  // Change the size of the array.
  void
  resize(int newsize);

  //////////
  // Changes size, makes exact if currently smaller...
  void
  setsize(int newsize);

  //////////
  // This will allocate an array which is equal to the number of
  // elements in the array and copy the values over.
  //
  // _size is not changed.
  //
  // newsize is an optional parameter which indicates the desired
  // size.  If newsize is greater than _size the new array will have
  // newsize elements in it (nalloc = newsize).  If newsize is less
  // than _size then _size elemets will be allocated (nalloc = _size).
  void
  trim(int newsize = 0);

  //////////
  // Initialize all elements of the array
  void
  initialize(const T& val);

  //////////
  // Get the array information
  T*
  get_objs();

};

template<class T>
Array1<T>::Array1(const Array1<T>& a)
{
  _size  = a._size;
  nalloc = _size;
  objs   = new T[_size];
  for (int i = 0; i < _size; i++) {
    objs[i] = a.objs[i];
  }
  nalloc            = _size;
  default_grow_size = a.default_grow_size;
}

template<class T>
Array1<T>&
Array1<T>::operator=(const Array1<T>& copy)
{
  if (&copy == this) {
    // Doing A=A, so don't do anything
    return (*this);
  }
  setsize(copy._size);
  for (int i = 0; i < _size; i++) {
    objs[i] = copy.objs[i];
  }
  default_grow_size = copy.default_grow_size;
  return (*this);
}

template<class T>
int
Array1<T>::operator==(const Array1<T>& a) const
{
  for (int i = 0; i < _size; i++) {
    if (objs[i] != a.objs[i]) {
      return false;
    }
  }

  return true;
}

template<class T>
int
Array1<T>::operator!=(const Array1<T>& a) const
{
  for (int i = 0; i < _size; i++) {
    if (objs[i] != a.objs[i]) {
      return true;
    }
  }

  return false;
}

template<class T>
Array1<T>::Array1(int size, int gs, int asize)
{
  ASSERT(size >= 0);
  default_grow_size = gs;
  _size             = size;
  if (size) {
    if (asize <= size) {
      objs   = new T[size];
      nalloc = _size;
    } else {
      objs   = new T[asize];
      nalloc = asize;
    }
  } else {
    if (asize > 0) {
      objs   = new T[asize];
      nalloc = asize;
    } else {
      objs   = 0;
      nalloc = 0;
    }
  }
}

template<class T>
Array1<T>::~Array1()
{
  if (objs) {
    delete[] objs;
  }
}

template<class T>
void
Array1<T>::grow(int count, int grow_size)
{
  int newsize = _size + count;
  if (newsize > nalloc) {
    // Reallocate...
    int gs1      = newsize >> 2;
    int gs       = gs1 > grow_size ? gs1 : grow_size;
    int newalloc = newsize + gs;
    T* newobjs   = new T[newalloc];
    if (objs) {
      for (int i = 0; i < _size; i++) {
        newobjs[i] = objs[i];
      }
      delete[] objs;
    }
    objs   = newobjs;
    nalloc = newalloc;
  }
  _size = newsize;
}

template<class T>
void
Array1<T>::reserve(int n)
{
  if (n > nalloc) {
    // Reallocate...
    T* newobjs = new T[n];
    if (objs) {
      for (int i = 0; i < _size; i++) {
        newobjs[i] = objs[i];
      }
      delete[] objs;
    }
    objs   = newobjs;
    nalloc = n;
  }
}

template<class T>
void
Array1<T>::add(const T& obj)
{
  grow(1, default_grow_size);
  objs[_size - 1] = obj;
}

template<class T>
void
Array1<T>::insert(int idx, const T& obj)
{
  grow(1, default_grow_size);
  for (int i = _size - 1; i > idx; i--) {
    objs[i] = objs[i - 1];
  }
  objs[idx] = obj;
}

template<class T>
void
Array1<T>::remove(int idx)
{
  _size--;
  for (int i = idx; i < _size; i++) {
    objs[i] = objs[i + 1];
  }
}

template<class T>
void
Array1<T>::remove_all()
{
  _size = 0;
}

template<class T>
void
Array1<T>::resize(int newsize)
{
  if (newsize > _size) {
    grow(newsize - _size);
  } else {
    _size = newsize;
  }
}

template<class T>
void
Array1<T>::setsize(int newsize)
{
  if (newsize > nalloc) {        // have to reallocate...
    T* newobjs = new T[newsize]; // make it exact!
    if (objs) {
      for (int i = 0; i < _size; i++) {
        newobjs[i] = objs[i];
      }
      delete[] objs;
    }
    objs   = newobjs;
    nalloc = newsize;
  }
  _size = newsize;
}

template<class T>
void
Array1<T>::trim(int newsize)
{
  if (newsize < _size || newsize <= 0) {
    newsize = _size;
  }
  if (newsize == nalloc) {
    // We already have the correct number allocated
    return;
  }
  T* newobjs = new T[newsize];
  if (objs) {
    // Copy the data
    for (int i = 0; i < _size; i++) {
      newobjs[i] = objs[i];
    }
    // Delete the old bit of memory
    delete[] objs;
  }
  objs   = newobjs;
  nalloc = newsize;
}

template<class T>
void
Array1<T>::initialize(const T& val)
{
  for (int i = 0; i < _size; i++) {
    objs[i] = val;
  }
}

template<class T>
T*
Array1<T>::get_objs()
{
  return objs;
}

} // End namespace Uintah

#endif /* __CORE_CONTAINERS_ARRAY1_H__ */
