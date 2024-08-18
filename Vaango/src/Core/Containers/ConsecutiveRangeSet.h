/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2015-2024 Biswajit Banerjee
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
#include <utility>
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
  ConsecutiveRangeSetException(std::string msg, const char* file, int line)
    : Exception()
    , msg_(std::move(msg))
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
    : Exception()
    , msg_(copy.msg_)
  {
  }

  auto
  operator=(const ConsecutiveRangeSetException& copy)
    -> ConsecutiveRangeSetException
  {
    msg_ = copy.msg_;
    return *this;
  }

  [[nodiscard]] const char*
  message() const override
  {
    return msg_.c_str();
  }

  [[nodiscard]] const char*
  type() const override
  {
    return "ConsecutiveRangeSetException";
  }

private:
  std::string msg_;
};

class ConsecutiveRangeSet
{
  friend auto
  operator<<(std::ostream& out, const ConsecutiveRangeSet& set)
    -> std::ostream&;

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

    iterator(const iterator& it2) = default;

    auto
    operator=(const iterator& it2) -> iterator& = default;

    inline auto
    operator*() noexcept(false) -> int;

    auto
    operator==(const iterator& it2) const -> bool
    {
      return range_ == it2.range_ && offset_ == it2.offset_;
    }

    auto
    operator!=(const iterator& it2) const -> bool
    {
      return !(*this == it2);
    }

    auto
    operator++() -> iterator&;

    inline auto
    operator++(int) -> iterator;

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

      = default;

    auto
    operator=(const Range& r2) -> Range& = default;

    auto
    operator==(const Range& r2) const -> bool
    {
      return low_ == r2.low_ && extent_ == r2.extent_;
    }

    auto
    operator!=(const Range& r2) const -> bool
    {
      return low_ != r2.low_ || extent_ != r2.extent_;
    }

    auto
    operator<(const Range& r2) const -> bool
    {
      return low_ < r2.low_;
    }

    inline void
    display(std::ostream& out) const;

    [[nodiscard]] auto
    high() const -> int
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

  ConsecutiveRangeSet(const ConsecutiveRangeSet& set2) = default;

  ~ConsecutiveRangeSet() = default;

  auto
  operator=(const ConsecutiveRangeSet& set2) -> ConsecutiveRangeSet& = default;

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

  auto
  operator==(const ConsecutiveRangeSet& set2) const -> bool;

  auto
  operator!=(const ConsecutiveRangeSet& set2) const -> bool
  {
    return !(*this == set2);
  }

  // obtain the intersection of two sets
  [[nodiscard]] auto
  intersected(const ConsecutiveRangeSet& set2) const -> ConsecutiveRangeSet;

  // obtain the union of two sets
  [[nodiscard]] auto
  unioned(const ConsecutiveRangeSet& set2) const -> ConsecutiveRangeSet;

  // Could implement binary search on the range set, but this
  // wasn't needed so I didn't do it.  Perhaps in the future if
  // needed. -- Wayne
  // I needed it, so I implemented it.  -- Bryan
  auto
  find(int n) -> iterator;

  [[nodiscard]] inline auto
  begin() const -> iterator
  {
    return iterator(this, 0, 0);
  }

  [[nodiscard]] inline auto
  end() const -> iterator
  {
    return iterator(this, (int)rangeSet_.size(), 0);
  }

  [[nodiscard]] auto
  size() const -> unsigned long
  {
    return size_;
  }

  [[nodiscard]] auto
  toString() const -> std::string;

  // return a space separated list of integers
  [[nodiscard]] auto
  expandedString() const -> std::string;

  // used for debugging
  auto
  getNumRanges() -> int
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

inline auto
ConsecutiveRangeSet::iterator::operator*() noexcept(false) -> int
{
  CHECKARRAYBOUNDS(range_, 0, (long)set_->rangeSet_.size());
  return set_->rangeSet_[range_].low_ + offset_;
}

inline auto
ConsecutiveRangeSet::iterator::operator++(int) -> ConsecutiveRangeSet::iterator
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
