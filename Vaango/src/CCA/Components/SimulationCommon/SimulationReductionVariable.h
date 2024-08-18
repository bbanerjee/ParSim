/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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

#ifndef __CCA_COMPONENTS_SIMULATIONCOMMON_SIMULATIONREDUCTIONVAR_H__
#define __CCA_COMPONENTS_SIMULATIONCOMMON_SIMULATIONREDUCTIONVAR_H__

#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <CCA/Ports/DataWarehouse.h>

namespace Uintah {

//class DataWarehouse;

class SimulationReductionVariable {
 public:
  SimulationReductionVariable(std::string name,
                              const TypeDescription* varType,
                              bool varActive = false) {
    // Construct the label.
    VarLabel* nonconstVar = VarLabel::create(name, varType);
    nonconstVar->isReductionTask(false);
    d_label = nonconstVar;

    d_active = varActive;

    setBenignValue();
    reset();
  };

  virtual ~SimulationReductionVariable() { VarLabel::destroy(d_label); }

  void
  setBenignValue() {
    d_bool_var.setBenignValue();
    d_min_var.setBenignValue();
    d_max_var.setBenignValue();
  }

  void
  reset() {
    d_reduction = false;

    d_count           = 0;
    d_overriddenValue = false;
  }

  // The set value call should be used before the reduction.
  void
  setValue(DataWarehouse* new_dw, bool val) {
    // Idiot proofing - If the reduction has occured do an override.
    if (d_reduction)
      overrideValue(new_dw, val);
    else {
      d_bool_value    = val;
      d_overrideValue = true;
    }
  }

  void
  setValue(DataWarehouse* new_dw, double val) {
    // Idiot proofing - If the reduction has occured do an override.
    if (d_reduction)
      overrideValue(new_dw, val);
    else {
      d_double_value  = val;
      d_overrideValue = true;
    }
  }

  // The override value call should be used after the reduction.
  void
  overrideValue(DataWarehouse* new_dw, bool val) {
    if (d_reduction) {
      d_overriddenValue = true;

      new_dw->override(bool_or_vartype(val), d_label);

      // Get the reduced value.
      new_dw->get(d_bool_var, d_label);
    }
    // Idiot proofing - If the reduction has not occured do a set.
    else
      setValue(new_dw, val);
  }

  void
  overrideValue(DataWarehouse* new_dw, double val) {
    if (d_reduction) {
      d_overriddenValue = true;

      if (d_label->typeDescription() == min_vartype::getTypeDescription()) {
        new_dw->override(min_vartype(val), d_label);
        // Get the reduced value.
        new_dw->get(d_max_var, d_label);
      } else if (d_label->typeDescription() ==
                 max_vartype::getTypeDescription()) {
        new_dw->override(max_vartype(val), d_label);
        // Get the reduced value.
        new_dw->get(d_max_var, d_label);
      }
    }
    // Idiot proofing - If the reduction has not occured do a set.
    else
      setValue(new_dw, val);
  }

  void
  reduce(DataWarehouse* new_dw) {
    Patch* patch = nullptr;

    d_bool_var.setBenignValue();
    d_min_var.setBenignValue();
    d_max_var.setBenignValue();

    // Reduce only if active.
    if (d_active) {
      if (d_label->typeDescription() == bool_or_vartype::getTypeDescription()) {
        // If the user gave a value use it.
        if (d_overrideValue) {
          new_dw->put(bool_or_vartype(d_bool_value), d_label);
          d_overrideValue   = false;
          d_overriddenValue = true;
        }
        // If the value does not exists put a benign value into the
        // warehouse.
        else if (!new_dw->exists(d_label, -1, patch)) {
          new_dw->put(d_bool_var, d_label);
        }

        // Only reduce if on more than one rank
        if (Parallel::getMPISize() > 1) {
          new_dw->reduceMPI(d_label, 0, 0, -1);
        }

        // Get the reduced value.
        new_dw->get(d_bool_var, d_label);
      } else if (d_label->typeDescription() ==
                 min_vartype::getTypeDescription()) {
        // If the user gave a value use it.
        if (d_overrideValue) {
          new_dw->put(min_vartype(d_double_value), d_label);
          d_overrideValue   = false;
          d_overriddenValue = true;
        }
        // If the value does not exists put a benign value into the
        // warehouse.
        else if (!new_dw->exists(d_label, -1, patch)) {
          new_dw->put(d_min_var, d_label);
        }

        // Only reduce if on more than one rank
        if (Parallel::getMPISize() > 1) {
          new_dw->reduceMPI(d_label, 0, 0, -1);
        }

        // Get the reduced value.
        new_dw->get(d_min_var, d_label);
      } else if (d_label->typeDescription() ==
                 max_vartype::getTypeDescription()) {
        // If the user gave a value use it.
        if (d_overrideValue) {
          new_dw->put(max_vartype(d_double_value), d_label);
          d_overrideValue = false;
        }
        // If the value does not exists put a benign value into the
        // warehouse.
        else if (!new_dw->exists(d_label, -1, patch)) {
          new_dw->put(d_max_var, d_label);
        }

        // Only reduce if on more than one rank
        if (Parallel::getMPISize() > 1) {
          new_dw->reduceMPI(d_label, 0, 0, -1);
        }

        // Get the reduced value.
        new_dw->get(d_max_var, d_label);
      }
    }

    d_reduction = true;
  }

  void
  setActive(bool val) {
    d_active = val;
  }
  bool
  getActive() const {
    return d_active;
  }

  const VarLabel*
  getLabel() const {
    return d_label;
  }
  const std::string
  getName() const {
    return d_label->getName();
  }

  double
  getValue() const {
    if (d_active) {
      if (d_label->typeDescription() == bool_or_vartype::getTypeDescription())
        return double(d_bool_var);
      else if (d_label->typeDescription() == min_vartype::getTypeDescription())
        return d_min_var;
      else if (d_label->typeDescription() == max_vartype::getTypeDescription())
        return d_max_var;
      else
        return 0;
    } else
      return 0;
  }

  bool
  isBenignValue() const {
    if (d_active) {
      if (d_label->typeDescription() == bool_or_vartype::getTypeDescription())
        return d_bool_var.isBenignValue();
      else if (d_label->typeDescription() == min_vartype::getTypeDescription())
        return d_min_var.isBenignValue();
      else if (d_label->typeDescription() == max_vartype::getTypeDescription())
        return d_max_var.isBenignValue();
      else
        return true;
    } else
      return true;
  }

  unsigned int
  getCount() const {
    return d_count;
  }
  bool
  overridden() const {
    return d_overriddenValue;
  }

 private:
  bool d_active{false};

  const VarLabel* d_label{nullptr};

  // Flag to indicate if the reduction has occured.
  bool d_reduction{false};
  // Count the number of times the value may have been set before
  // being cleared.
  unsigned int d_count{0};

  // Because this class gets put into a map it can not be a
  // template. As such, there are multiple storage variables. The
  // user need to know which one to use. Which they should given the
  // type description.

  // Also these vars hold the value of the reduction throughout the
  // execution of the task and do not get updated until the
  // reduction occurs. This scheme makes it possible to use the
  // value while the current time step is calculating the next
  // value.
  bool_or_vartype d_bool_var;
  min_vartype d_min_var;
  max_vartype d_max_var;

  // If the user has direct access to the application they can set
  // the reduction value directly rather than setting it in the data
  // warehouse. Before the reduction this value will get put into the
  // data warehouse for them thus overidding the current value.
  bool d_bool_value{0};
  double d_double_value{0};
  bool d_overrideValue{false};

  bool d_overriddenValue{false};
};

}  // End namespace Uintah

#endif  //__CCA_COMPONENTS_SIMULATIONCOMMON_SIMULATIONREDUCTIONVAR_H__
