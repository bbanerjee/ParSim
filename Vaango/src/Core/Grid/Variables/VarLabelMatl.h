/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#ifndef __CORE_GRID_VARIABLES_VARLABELMATL_H__
#define __CORE_GRID_VARIABLES_VARLABELMATL_H__

#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/VarLabel.h>

namespace Uintah {

/**************************************
  struct
    VarLabelMatl

    VarLabel, Material, and Domain


  GENERAL INFORMATION

    VarLabelMatl.h

    Wayne Witzel
    Department of Computer Science
    University of Utah

    Center for the Simulation of Accidental Fires and Explosions (C-SAFE)

  DESCRIPTION
    Specifies a VarLabel on a specific patch and for a specific material
    with an operator< defined so this can be used as a key in a map.
  ****************************************/

template<class DomainType0, class DomainType1 = void>
struct VarLabelMatl
{
  VarLabelMatl(const VarLabel* label,
               int mat_index,
               const DomainType0* domain0,
               const DomainType1* domain1)
    : d_label{ label }
    , d_matl_index{ mat_index }
    , d_domain0{ domain0 }
    , d_domain1{ domain1 }
  {
  }

  VarLabelMatl(const VarLabelMatl<DomainType0, DomainType1>& copy)
    : d_label{ copy.d_label }
    , d_matl_index{ copy.d_matl_index }
    , d_domain0{ copy.d_domain0 }
    , d_domain1{ copy.d_domain1 }
  {
  }

  VarLabelMatl<DomainType0, DomainType1>&
  operator=(const VarLabelMatl<DomainType0, DomainType1>& copy)
  {
    d_label      = copy.d_label;
    d_matl_index = copy.d_matl_index;
    d_domain0    = copy.d_domain0;
    d_domain1    = copy.d_domain1;
    return *this;
  }

  bool
  operator<(const VarLabelMatl<DomainType0, DomainType1>& other) const
  {
    if (d_label->equals(other.d_label)) {
      if (d_matl_index == other.d_matl_index) {
        if (d_domain0 == other.d_domain0) {
          return d_domain1 < other.d_domain1;
        } else {
          return d_domain0 < other.d_domain0;
        }
      } else {
        return d_matl_index < other.d_matl_index;
      }
    } else {
      VarLabel::Compare comp;
      return comp(d_label, other.d_label);
    }
  }

  bool
  operator==(const VarLabelMatl<DomainType0, DomainType1>& other) const
  {
    return ((d_label->equals(other.d_label)) &&
            (d_matl_index == other.d_matl_index) &&
            (d_domain0 == other.d_domain0) && (d_domain1 == other.d_domain1));
  }

  const VarLabel* d_label{ nullptr };
  int d_matl_index{};
  const DomainType0* d_domain0{ nullptr };
  const DomainType1* d_domain1{ nullptr };
};

template<class DomainType>
struct VarLabelMatl<DomainType, void>
{

  VarLabelMatl(const VarLabel* label, int matlIndex, const DomainType* domain)
    : m_label{ label }
    , m_matl_index{ matlIndex }
    , m_domain{ domain }
  {
  }

  VarLabelMatl(const VarLabelMatl<DomainType>& copy)
    : m_label(copy.m_label)
    , m_matl_index(copy.m_matl_index)
    , m_domain(copy.m_domain)
  {
  }

  VarLabelMatl<DomainType>&
  operator=(const VarLabelMatl<DomainType>& copy)
  {
    m_label      = copy.m_label;
    m_matl_index = copy.m_matl_index;
    m_domain     = copy.m_domain;
    return *this;
  }

  bool
  operator<(const VarLabelMatl<DomainType>& other) const
  {
    if (m_label->equals(other.m_label)) {
      if (m_matl_index == other.m_matl_index) {
        return m_domain < other.m_domain;
      } else {
        return m_matl_index < other.m_matl_index;
      }
    } else {
      VarLabel::Compare comp;
      return comp(m_label, other.m_label);
    }
  }

  bool
  operator==(const VarLabelMatl<DomainType>& other) const
  {
    return ((m_label->equals(other.m_label)) &&
            (m_matl_index == other.m_matl_index) &&
            (m_domain == other.m_domain));
  }

  const VarLabel* m_label{ nullptr };
  int m_matl_index{};
  const DomainType* m_domain{ nullptr };
};

} // End namespace Uintah

#endif //__CORE_GRID_VARIABLES_VARLABELMATL_H__
