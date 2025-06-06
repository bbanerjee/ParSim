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

#ifndef __CORE_GRID_VARIABLES_COMPUTESET_H__
#define __CORE_GRID_VARIABLES_COMPUTESET_H__

#include <Core/Exceptions/InternalError.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/Assert.h>
#include <Core/Util/RefCounted.h>
#include <Core/Util/constHandle.h>

#include <algorithm>
#include <sstream>
#include <vector>
#include <memory>

namespace Uintah {

/**************************************
  CLASS
  ComputeSet

  Provides similar functionality to std::set.  An exception is that
  a ComputeSet stores data in groups of ComputeSubsets.

  A ComputeSubset is much more similar to a std::set and provides
  functionality necessary for ComputeSet.

  GENERAL INFORMATION

  ComputeSet.h

  Steven G. Parker
  Department of Computer Science
  University of Utah

  Center for the Simulation of Accidental Fires and Explosions (C-SAFE)

 ****************************************/

class Patch;
template <class T>
class ComputeSubset;
template <class T>
class ComputeSet;

typedef ComputeSet<const Patch*> PatchSet;
typedef ComputeSet<int> MaterialSet;
typedef ComputeSubset<const Patch*> PatchSubset;
typedef ComputeSubset<int> MaterialSubset;

using MaterialSetSP = std::shared_ptr<MaterialSet>;
using MaterialSubsetSP = std::shared_ptr<MaterialSubset>;

template <class T>
class ComputeSubset : public RefCounted
{
public:
  ComputeSubset(int size = 0)
    : d_items(size)
  {
  }
  ComputeSubset(const std::vector<T>& items)
    : d_items(items)
  {
  }

  template <class InputIterator>
  ComputeSubset(InputIterator begin, InputIterator end)
    : d_items(begin, end)
  {
  }

  auto begin() { return d_items.begin(); }
  auto end() { return d_items.end(); }
  const auto begin() const { return d_items.begin(); }
  const auto end() const { return d_items.end(); }

  int size() const { return (int)d_items.size(); }

  T& operator[](int i) { return d_items[i]; }
  const T& operator[](int i) const { return d_items[i]; }
  const T& get(int i) const { return d_items[i]; }
  void add(const T& i) { d_items.push_back(i); }

  void addSubset(const ComputeSubset<T>* subset);

  bool empty() const { return d_items.size() == 0; }
  bool contains(T elem) const
  {
    for (unsigned int i = 0; i < d_items.size(); i++) {
      if (d_items[i] == elem) {
        return true;
      }
    }
    return false;
  }

  void sort();
  bool is_sorted() const;

  const std::vector<T>& getVector() const { return d_items; }

  bool equals(const ComputeSubset<T>* s2) const
  {
    // check that the sets are equvalent
    if (d_items.size() != s2->d_items.size()) {
      return false;
    }

    for (unsigned int i = 0; i < d_items.size(); i++) {
      if (d_items[i] != s2->d_items[i]) {
        return false;
      }
    }

    return true;
  };

  constHandle<ComputeSubset<T>> intersection(
    constHandle<ComputeSubset<T>> s2) const
  {
    return intersection(this, s2);
  }

  // May return the same Handle as one that came in.
  static constHandle<ComputeSubset<T>> intersection(
    const constHandle<ComputeSubset<T>>& s1,
    const constHandle<ComputeSubset<T>>& s2);

  // May pass back Handles to same sets that came in.
  static void intersectionAndDifferences(
    const constHandle<ComputeSubset<T>>& A,
    const constHandle<ComputeSubset<T>>& B,
    constHandle<ComputeSubset<T>>& intersection,
    constHandle<ComputeSubset<T>>& AminusB,
    constHandle<ComputeSubset<T>>& BminusA)
  {
    intersection =
      intersectionAndMaybeDifferences<true>(A, B, AminusB, BminusA);
  }

  static void difference(const constHandle<ComputeSubset<T>>& A,
                         const constHandle<ComputeSubset<T>>& B,
                         constHandle<ComputeSubset<T>>& diff)
  {
    diff = difference(A, B);
  }

  static bool overlaps(const ComputeSubset<T>* s1, const ComputeSubset<T>* s2);

  static bool compareElems(const T e1, const T e2);

private:
  // May pass back Handles to same sets that came in.
  template <bool passBackDifferences>
  static constHandle<ComputeSubset<T>> intersectionAndMaybeDifferences(
    const constHandle<ComputeSubset<T>>& s1,
    const constHandle<ComputeSubset<T>>& s2,
    constHandle<ComputeSubset<T>>& setDifference1,
    constHandle<ComputeSubset<T>>& setDifference2);

  static constHandle<ComputeSubset<T>> difference(
    const constHandle<ComputeSubset<T>>& A,
    const constHandle<ComputeSubset<T>>& B);

  std::vector<T> d_items;

  ComputeSubset(const ComputeSubset&) = delete;
  ComputeSubset& operator=(const ComputeSubset&) = delete;

}; // end class ComputeSubset

template <class T>
class ComputeSet : public RefCounted
{
public:
  ComputeSet();
  ~ComputeSet();

  // adds all unique elements of vector in one subset
  void addAll_unique(const std::vector<T>&);

  // adds all elements of vector in one subset
  void addAll(const std::vector<T>&);

  // adds each element of vector as a separate individual subset
  void addEach(const std::vector<T>&);

  // adds one element as a new subset
  void add(const T&);

  // add existing subset to the set
  void addSubset(ComputeSubset<T>* subset);

  void sortSubsets();

  int size() const { return (int)d_set.size(); }
  ComputeSubset<T>* getSubset(int idx) { return d_set[idx]; }

  /// Returns the vector of subsets managed by this set
  const std::vector<ComputeSubset<T>*>& getVector() const { return d_set; }

  const ComputeSubset<T>* getSubset(int idx) const { return d_set[idx]; }
  const ComputeSubset<T>* getUnion() const;
  void createEmptySubsets(int size);
  int totalsize() const;

private:
  std::vector<ComputeSubset<T>*> d_set;
  mutable ComputeSubset<T>* d_subset;

  ComputeSet(const ComputeSet&) = delete;
  ComputeSet& operator=(const ComputeSet&) = delete;

  friend std::ostream& operator<<(std::ostream& out, const Uintah::PatchSet&);
  friend std::ostream& operator<<(std::ostream& out,
                                  const Uintah::MaterialSet&);
  friend std::ostream& operator<<(std::ostream& out,
                                  const Uintah::PatchSubset&);
  friend std::ostream& operator<<(std::ostream& out,
                                  const Uintah::MaterialSubset&);

}; // end class ComputeSet

std::ostream& operator<<(std::ostream& out, const Uintah::PatchSet&);
std::ostream& operator<<(std::ostream& out, const Uintah::MaterialSet&);
std::ostream& operator<<(std::ostream& out, const Uintah::PatchSubset&);
std::ostream& operator<<(std::ostream& out, const Uintah::MaterialSubset&);

template <class T>
ComputeSet<T>::ComputeSet()
{
  d_subset = nullptr;
}

template <class T>
ComputeSet<T>::~ComputeSet()
{
  for (unsigned int i = 0; i < d_set.size(); i++) {
    if (d_set[i]->removeReference())
      delete d_set[i];
  }
  if (d_subset && d_subset->removeReference()) {
    delete d_subset;
  }
}

template <class T>
void
ComputeSet<T>::addAll(const std::vector<T>& sub)
{
  ASSERT(!d_subset);
  ComputeSubset<T>* subset = scinew ComputeSubset<T>(sub);
  subset->sort();
  subset->addReference();
  d_set.push_back(subset);
}

template <class T>
void
ComputeSet<T>::addAll_unique(const std::vector<T>& sub)
{
  // Only insert unique entries into the set
  std::vector<T> sub_unique;
  sub_unique.push_back(sub[0]);

  for (unsigned int i = 1; i < sub.size(); i++) {
    if (find(sub_unique.begin(), sub_unique.end(), sub[i]) ==
        sub_unique.end()) {
      sub_unique.push_back(sub[i]);
    }
  }

  ASSERT(!d_subset);
  ComputeSubset<T>* subset = scinew ComputeSubset<T>(sub_unique);
  subset->sort();
  subset->addReference();
  d_set.push_back(subset);
}

template <class T>
void
ComputeSet<T>::addEach(const std::vector<T>& sub)
{
  ASSERT(!d_subset);
  for (unsigned int i = 0; i < sub.size(); i++) {
    ComputeSubset<T>* subset = scinew ComputeSubset<T>(1);
    subset->addReference();
    (*subset)[0] = sub[i];
    d_set.push_back(subset);
  }
}

template <class T>
void
ComputeSet<T>::add(const T& item)
{
  ASSERT(!d_subset);
  ComputeSubset<T>* subset = scinew ComputeSubset<T>(1);
  subset->addReference();
  (*subset)[0] = item;
  d_set.push_back(subset);
}

template <class T>
void
ComputeSet<T>::addSubset(ComputeSubset<T>* sub)
{
  ASSERT(!d_subset);
  ComputeSubset<T>* subset = sub; // hack to get around const issue
  d_set.push_back(subset);
}

template <class T>
void
ComputeSet<T>::sortSubsets()
{
  for (unsigned int i = 0; i < d_set.size(); i++) {
    ComputeSubset<T>* ss = d_set[i];
    ss->sort();
  }
}

template <class T>
void
ComputeSet<T>::createEmptySubsets(int n)
{
  ASSERT(!d_subset);
  for (int i = 0; i < n; i++) {
    ComputeSubset<T>* subset = scinew ComputeSubset<T>(0);
    subset->addReference();
    d_set.push_back(subset);
  }
}

template <class T>
const ComputeSubset<T>*
ComputeSet<T>::getUnion() const
{
  if (!d_subset) {
    d_subset = scinew ComputeSubset<T>;
    d_subset->addReference();
    for (unsigned int i = 0; i < d_set.size(); i++) {
      ComputeSubset<T>* ss = d_set[i];
      for (int j = 0; j < ss->size(); j++) {
        d_subset->add(ss->get(j));
      }
    }
    d_subset->sort();
  }
  return d_subset;
}

template <class T>
int
ComputeSet<T>::totalsize() const
{
  if (d_subset) {
    return d_subset->size();
  }
  int total = 0;
  for (unsigned int i = 0; i < d_set.size(); i++)
    total += d_set[i]->size();
  return total;
}

// Note: sort is specialized in ComputeSet_special for const Patch*'s
// to use Patch::Compare.
template <>
void ComputeSubset<const Patch*>::sort();

template <class T>
void
ComputeSubset<T>::sort()
{
  std::sort(d_items.begin(), d_items.end());
}

// specialized for patch in ComputeSet_special.cc
template <>
bool ComputeSubset<const Patch*>::compareElems(const Patch* e1,
                                               const Patch* e2);

template <class T>
bool
ComputeSubset<T>::compareElems(const T e1, const T e2)
{
  return e1 < e2;
}

template <class T>
void
ComputeSubset<T>::addSubset(const ComputeSubset<T>* ss)
{
  for (auto m = 0; m < ss->size(); m++) {
    add(ss->get(m));
  }
}

// compute the interesection between s1 and s2
template <class T>
constHandle<ComputeSubset<T>>
ComputeSubset<T>::intersection(const constHandle<ComputeSubset<T>>& s1,
                               const constHandle<ComputeSubset<T>>& s2)
{
  if (s1 == s2) {
    // for efficiency -- expedite when s1 and s2 point to the same thing
    return s1;
  }

  if (!s1) {
    return s2; // treat null as everything -- intersecting with it gives all
  }
  if (!s2) {
    return s1; // treat null as everything -- intersecting with it gives all
  }

  if (s1->size() == 0)
    return s1; // return an empty set
  if (s2->size() == 0)
    return s2; // return an empty set

  Handle<ComputeSubset<T>> intersection = scinew ComputeSubset<T>;

#if SCI_ASSERTION_LEVEL > 0
  if (!s1->is_sorted()) {
    SCI_THROW(InternalError("ComputeSubset s1 not sorted in "
                            "ComputeSubset<T>::intersectionAndMaybeDifference",
                            __FILE__, __LINE__));
  }
  if (!s2->is_sorted()) {
    SCI_THROW(InternalError("ComputeSubset s2 not sorted in "
                            "ComputeSubset<T>::intersectionAndMaybeDifference",
                            __FILE__, __LINE__));
  }

  T el2 = s2->get(0);
  for (int i = 1; i < s2->size(); i++) {
    T el = s2->get(i);
    if (!compareElems(el2, el)) {
      std::ostringstream msgstr;
      msgstr << "Set not sorted: " << el2 << ", " << el;
      SCI_THROW(InternalError(msgstr.str(), __FILE__, __LINE__));
    }
    el2 = el;
  }
#endif
  int i1 = 0;
  int i2 = 0;
  for (;;) {
    if (s1->get(i1) == s2->get(i2)) {
      intersection->add(s1->get(i1));
      i1++;
      i2++;
    } else if (compareElems(s1->get(i1), s2->get(i2))) {
      i1++;
    } else {
      i2++;
    }
    if (i1 == s1->size() || i2 == s2->size())
      break;
  }

  return intersection;
}

template <class T>
template <bool passBackDifferences>
constHandle<ComputeSubset<T>>
ComputeSubset<T>::intersectionAndMaybeDifferences(
  const constHandle<ComputeSubset<T>>& s1,
  const constHandle<ComputeSubset<T>>& s2,
  constHandle<ComputeSubset<T>>& setDifference1,
  constHandle<ComputeSubset<T>>& setDifference2)
{
  if (s1 == s2) {
    // for efficiency -- expedite when s1 and s2 point to the same thing
    setDifference1 = setDifference2 = scinew ComputeSubset<T>(0);
    return s1;
  }

  if (passBackDifferences) {
    setDifference1 = s1;
    setDifference2 = s2;
  }

  if (!s1) {
    setDifference2 = 0; // arbitrarily
    return s2; // treat null as everything -- intersecting with it gives all
  }
  if (!s2) {
    setDifference1 = 0; // arbitrarily
    return s1; // treat null as everything -- intersecting with it gives all
  }

  if (s1->size() == 0)
    return s1; // return an empty set
  if (s2->size() == 0)
    return s2; // return an empty set

  Handle<ComputeSubset<T>> intersection = scinew ComputeSubset<T>;
  Handle<ComputeSubset<T>> s1_minus_s2, s2_minus_s1;

  if (passBackDifferences) {
    setDifference1 = s1_minus_s2 = scinew ComputeSubset<T>;
    setDifference2 = s2_minus_s1 = scinew ComputeSubset<T>;
  }

#if SCI_ASSERTION_LEVEL > 0
  if (!s1->is_sorted()) {
    SCI_THROW(InternalError("ComputeSubset s1 not sorted in "
                            "ComputeSubset<T>::intersectionAndMaybeDifference",
                            __FILE__, __LINE__));
  }
  if (!s2->is_sorted()) {
    SCI_THROW(InternalError("ComputeSubset s2 not sorted in "
                            "ComputeSubset<T>::intersectionAndMaybeDifference",
                            __FILE__, __LINE__));
  }

  T el2 = s2->get(0);
  for (int i = 1; i < s2->size(); i++) {
    T el = s2->get(i);
    if (!compareElems(el2, el)) {
      std::ostringstream msgstr;
      msgstr << "Set not sorted: " << el2 << ", " << el;
      SCI_THROW(InternalError(msgstr.str(), __FILE__, __LINE__));
    }
    el2 = el;
  }
#endif
  int i1 = 0;
  int i2 = 0;
  for (;;) {
    if (s1->get(i1) == s2->get(i2)) {
      intersection->add(s1->get(i1));
      i1++;
      i2++;
    } else if (compareElems(s1->get(i1), s2->get(i2))) {
      if (passBackDifferences) {
        s1_minus_s2->add(s1->get(i1)); // alters setDifference1
      }
      i1++;
    } else {
      if (passBackDifferences) {
        s2_minus_s1->add(s2->get(i2)); // alters setDifference2
      }
      i2++;
    }
    if (i1 == s1->size() || i2 == s2->size())
      break;
  }

  if (passBackDifferences) {
    if (intersection->empty()) {
      // if the intersection is empty, then the set differences are
      // the same as the sets that came in.
      setDifference1 = s1;
      setDifference2 = s2;
    } else {
      // get the rest of whichever difference set wasn't finished (if any)
      for (; i1 != s1->size(); i1++) {
        s1_minus_s2->add(s1->get(i1)); // alters setDifference1
      }
      for (; i2 != s2->size(); i2++) {
        s2_minus_s1->add(s2->get(i2)); // alters setDifference2
      }
    }
  }

  if (intersection->size() == s1->size()) {
    return s1;
  } else if (intersection->size() == s2->size()) {
    return s2;
  } else {
    return intersection;
  }
}

template <class T>
constHandle<ComputeSubset<T>>
ComputeSubset<T>::difference(const constHandle<ComputeSubset<T>>& s1,
                             const constHandle<ComputeSubset<T>>& s2)
{
  if (s1 == s2 || !s1 || !s2) {
    // for efficiency -- expedite when s1 and s2 point to the same thing or are
    // null
    return Handle<ComputeSubset<T>>(scinew ComputeSubset<T>);
  }

  if (s1->size() == 0 || s1->size() == 0 || s2->size() == 0) {
    return s1; // return s1
  }

#if SCI_ASSERTION_LEVEL > 0
  if (!s1->is_sorted()) {
    SCI_THROW(InternalError("ComputeSubset s1 not sorted in "
                            "ComputeSubset<T>::intersectionAndMaybeDifference",
                            __FILE__, __LINE__));
  }
  if (!s2->is_sorted()) {
    SCI_THROW(InternalError("ComputeSubset s2 not sorted in "
                            "ComputeSubset<T>::intersectionAndMaybeDifference",
                            __FILE__, __LINE__));
  }

  T el2 = s2->get(0);
  for (int i = 1; i < s2->size(); i++) {
    T el = s2->get(i);
    if (!compareElems(el2, el)) {
      std::ostringstream msgstr;
      msgstr << "Set not sorted: " << el2 << ", " << el;
      SCI_THROW(InternalError(msgstr.str(), __FILE__, __LINE__));
    }
    el2 = el;
  }
#endif

  Handle<ComputeSubset<T>> diff = scinew ComputeSubset<T>;

  int i1 = 0;
  int i2 = 0;
  for (;;) {
    if (s1->get(i1) == s2->get(i2)) {
      // in both s1 and s2
      i1++;
      i2++;
    } else if (compareElems(s1->get(i1), s2->get(i2))) {
      // in s1 but not s2
      diff->add(s1->get(i1)); // alters setDifference1
      i1++;
    } else {
      // in s2 but not s1
      i2++;
    }
    if (i1 == s1->size() || i2 == s2->size())
      break;
  }

  // get the rest of s1 that wasn't finished (if any)
  for (; i1 != s1->size(); i1++) {
    diff->add(s1->get(i1)); // alters setDifference1
  }

  return diff;
}

template <class T>
bool
ComputeSubset<T>::overlaps(const ComputeSubset<T>* s1,
                           const ComputeSubset<T>* s2)
{
  if (s1 == s2) {
    return true;
  }
  if (s1->size() == 0 || s2->size() == 0) {
    return false;
  }
  if (!s1->is_sorted()) {
    SCI_THROW(
      InternalError("ComputeSubset s1 not sorted in ComputeSubset<T>::overlaps",
                    __FILE__, __LINE__));
  }
  if (!s2->is_sorted()) {
    SCI_THROW(
      InternalError("ComputeSubset s2 not sorted in ComputeSubset<T>::overlaps",
                    __FILE__, __LINE__));
  }
  int i1 = 0;
  int i2 = 0;
  for (;;) {
    if (s1->get(i1) == s2->get(i2)) {
      return true;
    } else if (compareElems(s1->get(i1), s2->get(i2))) {
      if (++i1 == s1->size())
        break;
    } else {
      if (++i2 == s2->size())
        break;
    }
  }
  return false;
}

template <class T>
bool
ComputeSubset<T>::is_sorted() const
{
  for (int i = 1; i < size(); ++i) {
    if (!compareElems(get(i - 1), get(i))) {
      return false;
    }
  }
  return true;
}

} // end namespace Uintah

#endif // _CORE_GRID_VARIABLES_COMPUTESET_H__
