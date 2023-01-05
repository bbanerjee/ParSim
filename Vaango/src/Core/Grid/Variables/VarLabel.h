/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2018-2023 Parresia Research Limited, NZ
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

#ifndef __CORE_GRID_VARIABLES_VARLABEL_H__
#define __CORE_GRID_VARIABLES_VARLABEL_H__

#include <Core/Geometry/IntVector.h>
#include <Core/Util/RefCounted.h>

#include <iosfwd>
#include <map>
#include <string>

namespace Uintah {

class TypeDescription;
class Patch;

class VarLabel : public RefCounted {
 public:
  enum class VarType { Normal, PositionVariable };

  // Ensure the uniqueness of VarLabel names (same name, same object).
  static VarLabel*
  create(const std::string& name,
         const TypeDescription* type_description,
         const IntVector& boundaryLayer = IntVector(0, 0, 0),
         VarType vartype                = VarType::Normal);

  static bool
  destroy(const VarLabel* label);

  inline const std::string&
  getName() const {
    return d_name;
  }

  std::string
  getFullName(int matlIndex, const Patch* patch) const;

  bool
  isPositionVariable() const {
    return d_var_type == VarType::PositionVariable;
  }

  const TypeDescription*
  typeDescription() const {
    return d_td;
  }

  IntVector
  getBoundaryLayer() const {
    return d_boundary_layer;
  }

  // void allowMultipleComputes();

  // bool allowsMultipleComputes() const { return d_allow_multiple_computes; }

  void
  isReductionTask(bool input);

  bool
  isReductionTask() const {
    return d_is_reduction_task;
  }

  static VarLabel*
  find(const std::string& name);

  static VarLabel*
  find(const std::string& name, const std::string& message);

  static VarLabel*
  particlePositionLabel();

  static void
  setParticlePositionName(const std::string& pPosName) {
    s_particle_position_name = pPosName;
  }

  static std::string&
  getParticlePositionName() {
    return s_particle_position_name;
  }

  class Compare {
   public:
    inline bool
    operator()(const VarLabel* v1, const VarLabel* v2) const {
      // because of uniqueness, we can use pointer comparisons
      // return v1 < v2;
      // No we cannot, because we need the order to be the same on different
      // processes
      if (v1 == v2) {
        return false;
      }
      return v1->getName() < v2->getName();
    }
  };

  bool
  equals(const VarLabel* v2) const {
    // because of uniqueness, we can use pointer comparisons
    return this == v2;
    /* old way
       if(this == v2)
       return true;
       return getName() == v2->getName();
    */
  }

  void
  setCompressionMode(std::string compressionMode) {
    d_compression_mode = compressionMode;
  }

  const std::string&
  getCompressionMode() const {
    return (d_compression_mode == "default") ? s_default_compression_mode
                                             : d_compression_mode;
  }

  static void
  setDefaultCompressionMode(const std::string& compressionMode) {
    s_default_compression_mode = compressionMode;
  }

  static void
  printAll();  // for debugging

  friend std::ostream&
  operator<<(std::ostream& out, const VarLabel& vl);

 private:
  // You must use VarLabel::create.
  VarLabel(const std::string& label,
           const TypeDescription* type,
           const IntVector& boundaryLayer,
           VarType vartype);

  // You must use destroy.
  ~VarLabel(){};

 private:
  std::string d_name{""};
  const TypeDescription* d_td{nullptr};
  IntVector d_boundary_layer{IntVector(0, 0, 0)};
  VarType d_var_type{VarType::Normal};

  mutable std::string d_compression_mode{"default"};
  static std::string s_default_compression_mode;
  static std::string s_particle_position_name;

  // Allow a variable of this label to be computed multiple times in a TaskGraph
  // without complaining.
  // bool d_allow_multiple_computes{ false };
  bool d_is_reduction_task{true};

  // eliminate copy, assignment and move
  VarLabel(const VarLabel&) = delete;
  VarLabel&
  operator=(const VarLabel&) = delete;
  VarLabel(VarLabel&&)       = delete;
  VarLabel&
  operator=(VarLabel&&) = delete;

  // Static member to keep track of all labels created to prevent
  // duplicates.
  static std::map<std::string, VarLabel*> g_all_labels;
};

}  // End namespace Uintah

#endif  // __CORE_GRID_VARIABLES_VARLABEL_H__
