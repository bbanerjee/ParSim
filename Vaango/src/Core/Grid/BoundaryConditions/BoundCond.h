/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2014-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

#ifndef __CORE_GRID_BOUNDARYCONDITIONS_BoundCond_H__
#define __CORE_GRID_BOUNDARYCONDITIONS_BoundCond_H__

#include <Core/Grid/BoundaryConditions/BoundCondBase.h>

#include <Core/Geometry/Vector.h>
#include <Core/Malloc/Allocator.h>

#include <string>

namespace Uintah {

class NoValue
{
public:
  NoValue()  = default;
  ~NoValue() = default;
};

} // namespace Uintah

namespace Uintah {

template<class T>
class BoundCond : public BoundCondBase
{
public:
  using BoundCondP = std::shared_ptr<BoundCond<T>>;

  BoundCond() = default;

  BoundCond(const std::string& var_name,
            const std::string& type,
            const T value,
            const std::string& face_label,
            const BoundCondBase::BoundCondValueTypeEnum val_type)
    : d_value{ value }
  {
    d_variable   = var_name;
    d_type       = type;
    d_face_label = face_label;
    d_value_type = val_type;
  }

  ~BoundCond() override = default;

  auto
  clone() -> BoundCondP
  {
    return BoundCondP(cloneImpl());
  };

  auto
  getValue() const -> T
  {
    return d_value;
  };

  [[nodiscard]] auto
  getType() const -> const std::string
  {
    return d_type;
  };

protected:
  std::shared_ptr<BoundCondBase>
  cloneImpl() override
  {
    return std::make_shared<BoundCond>(*this);
  }

protected:
  T d_value;
};

} // namespace Uintah

namespace Uintah {

template<>
class BoundCond<NoValue> : public BoundCondBase
{

public:
  using BoundCondP = std::shared_ptr<BoundCond<NoValue>>;

public:
  BoundCond(string var_name, string type)
  {
    d_value    = NoValue();
    d_variable = var_name;
    d_type     = type;
  }

  BoundCond(string var_name)
  {
    d_value    = NoValue();
    d_variable = var_name;
  }

  auto
  clone() 
  {
    return cloneImpl();
  }

protected:
  std::shared_ptr<BoundCondBase> 
  cloneImpl() override
  {
    return std::make_shared<BoundCond>(*this);
  }

protected:
  NoValue d_value;
};

} // End namespace Uintah

#endif //__CORE_GRID_BOUNDARYCONDITIONS_BoundCond_H__
