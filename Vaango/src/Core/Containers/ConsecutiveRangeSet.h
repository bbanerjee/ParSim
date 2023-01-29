/*
 * The MIT License
 *
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

/*
 *  ConsecutiveRangeSet.h
 *
 *  Written by:
 *   Wayne Witzel
 *   Department of Computer Science
 *   University of Utah
 *   Nov. 2000
 *
 *
 */

#ifndef __CORE_CONTAINERS_CONSECUTIVE_RANGESET_H__
#define __CORE_CONTAINERS_CONSECUTIVE_RANGESET_H__

#include <Core/Exceptions/Exception.h>
#include <Core/Util/Assert.h>

#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

namespace Uintah {

/**************************************
CLASS
   ConsecutiveRangeSet

KEYWORDS
   Interval, integers

DESCRIPTION

   Represents a set of integers that are
   stored efficiently if grouped together
   in consecutive ranges.

   Written by:
    Wayne Witzel
    Department of Computer Science
    University of Utah
    Nov. 2000
****************************************/

class ConsecutiveRangeSetException : public Exception
{
public:
  ConsecutiveRangeSetException(const std::string& msg,
                               const char* file,
                               int line)
    : Exception(), msg_(msg)
  {

    std::ostringstream s;
    s << "A ConsecutiveRangeSetException exception was thrown\n"
      << file << ":" << line << "\n"
      << msg_;
    msg_ = (char*)(s.str().c_str());

#ifdef EXCEPTIONS_CRASH
    std::cout << msg_ << "\n";
#endif
  }

  ConsecutiveRangeSetException(const ConsecutiveRangeSetException& copy)
    : Exception(), msg_(copy.msg_)
  {
  }

  ConsecutiveRangeSetException
  operator=(const ConsecutiveRangeSetException& copy)
  {
    msg_ = copy.msg_;
    return *this;
  }

  virtual const char*
  message() const
  {
    return msg_.c_str();
  }

  virtual const char*
  type() const
  {
    return "ConsecutiveRangeSetException";
  }

private:
  std::string msg_;
};

class ConsecutiveRangeSet
{
  friend std::ostream&
  operator<<(std::ostream& out, const ConsecutiveRangeSet& set);

public:
  class iterator
  {
  public:
    iterator(const ConsecutiveRangeSet* set, int range, int offset)
      : set_{ set }
      , range_{ range }
      , offset_{ offset }
    {
    }

    iterator(const iterator& it2)
      : set_{ it2.set_ }
      , range_{ it2.range_ }
      , offset_{ it2.offset_ }
    {
    }

    iterator&
    operator=(const iterator& it2)
    {
      set_    = it2.set_;
      range_  = it2.range_;
      offset_ = it2.offset_;
      return *this;
    }

    inline int
    operator*() noexcept(false);

    bool
    operator==(const iterator& it2) const
    {
      return range_ == it2.range_ && offset_ == it2.offset_;
    }

    bool
    operator!=(const iterator& it2) const
    {
      return !(*this == it2);
    }

    iterator&
    operator++();

    inline iterator
    operator++(int);

  private:
    const ConsecutiveRangeSet* set_{ nullptr };
    int range_{ 0 };
    int offset_{ 0 };
  };

  // represents range: [low, low+extent]
  struct Range
  {
    Range(int low, int high);

    Range(const Range& r2)
      : low_{ r2.low_ }
      , extent_{ r2.extent_ }
    {
    }

    Range&
    operator=(const Range& r2)
    {
      low_    = r2.low_;
      extent_ = r2.extent_;
      return *this;
    }

    bool
    operator==(const Range& r2) const
    {
      return low_ == r2.low_ && extent_ == r2.extent_;
    }

    bool
    operator!=(const Range& r2) const
    {
      return low_ != r2.low_ || extent_ != r2.extent_;
    }

    bool
    operator<(const Range& r2) const
    {
      return low_ < r2.low_;
    }

    inline void
    display(std::ostream& out) const;

    int
    high() const
    {
      return (int)(low_ + extent_);
    }

    int low_{ 0 };
    unsigned long extent_{ 0 };
  };

public:
  ConsecutiveRangeSet(std::list<int>& set);

  ConsecutiveRangeSet(int low, int high); // single consecutive range

  // empty set
  ConsecutiveRangeSet() = default;

  // initialize a range set with a string formatted like: "1, 2-8, 10, 15-30"
  ConsecutiveRangeSet(const std::string& setstr) noexcept(false);

  ConsecutiveRangeSet(const ConsecutiveRangeSet& set2)
    : rangeSet_{ set2.rangeSet_ }
    , size_{ set2.size_ }
  {
  }

  ~ConsecutiveRangeSet() = default;

  ConsecutiveRangeSet&
  operator=(const ConsecutiveRangeSet& set2)
  {
    rangeSet_ = set2.rangeSet_;
    size_     = set2.size_;
    return *this;
  }

  // Add to the range set, asserting that value is greater or equal
  // to anything already in the set (if it is equal to something already
  // in teh set then the value is simply discarded).
  void
  addInOrder(int value) noexcept(false);

  template<class AnyIterator>
  void
  addInOrder(const AnyIterator& begin, const AnyIterator& end)
  {
    for (AnyIterator it = begin; it != end; ++it) {
      addInOrder(*it);
    }
  }

  bool
  operator==(const ConsecutiveRangeSet& set2) const;

  bool
  operator!=(const ConsecutiveRangeSet& set2) const
  {
    return !(*this == set2);
  }

  // obtain the intersection of two sets
  ConsecutiveRangeSet
  intersected(const ConsecutiveRangeSet& set2) const;

  // obtain the union of two sets
  ConsecutiveRangeSet
  unioned(const ConsecutiveRangeSet& set2) const;

  // Could implement binary search on the range set, but this
  // wasn't needed so I didn't do it.  Perhaps in the future if
  // needed. -- Wayne
  // I needed it, so I implemented it.  -- Bryan
  iterator
  find(int n);

  inline iterator
  begin() const
  {
    return iterator(this, 0, 0);
  }

  inline iterator
  end() const
  {
    return iterator(this, (int)rangeSet_.size(), 0);
  }

  unsigned long
  size() const
  {
    return size_;
  }

  std::string
  toString() const;

  // return a space separated list of integers
  std::string
  expandedString() const;

  // used for debugging
  int
  getNumRanges()
  {
    return (int)rangeSet_.size();
  }

  static const ConsecutiveRangeSet empty;
  static const ConsecutiveRangeSet all;
  friend class ConsecutiveRangeSet::iterator;

private:
  template<class InputIterator>
  ConsecutiveRangeSet(InputIterator begin, InputIterator end)
    : rangeSet_{ begin, end }
  {
    setSize();
  }
  void
  setSize();

  std::vector<Range> rangeSet_{};
  unsigned long size_{ 0 }; // sum of range (extent+1)'s
};

inline int
ConsecutiveRangeSet::iterator::operator*() noexcept(false)
{
  CHECKARRAYBOUNDS(range_, 0, (long)set_->rangeSet_.size());
  return set_->rangeSet_[range_].low_ + offset_;
}

inline ConsecutiveRangeSet::iterator
ConsecutiveRangeSet::iterator::operator++(int)
{
  iterator oldit(*this);
  ++(*this);
  return oldit;
}

inline void
ConsecutiveRangeSet::Range::display(std::ostream& out) const
{
  if (extent_ == 0) {
    out << low_;
  } else {
    out << low_ << " - " << high();
  }
}

} // End namespace Uintah

#endif //__CORE_CONTAINERS_CONSECUTIVE_RANGESET_H__
