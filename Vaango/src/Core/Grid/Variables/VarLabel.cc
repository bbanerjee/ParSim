/*
 * The MIT License
 *
 * Copyright (c) 1997-2022 The University of Utah
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

#include <Core/Grid/Variables/VarLabel.h>

#include <Core/Exceptions/InternalError.h>
#include <Core/Grid/Patch.h>
#include <Core/Parallel/MasterLock.h>
#include <Core/Util/DebugStream.h>

//#include <Core/Util/DebugOutput.h>

#include <iostream>
#include <map>
#include <memory>
#include <sstream>

namespace Uintah {

namespace {
MasterLock g_label_mutex{};
//DebugOutput g_varlabel_dbg("VarLabel", "VarLabel", 
//                           "Report Varlabel creation and deletion", false);
}

static DebugStream dbg("VarLabel", false);

//______________________________________________________________________
// Initialize class static variables:

std::string VarLabel::s_particle_position_name = "p.x";
std::string VarLabel::s_default_compression_mode = "none";

std::map<std::string, VarLabel*> VarLabel::g_all_labels;

//______________________________________________________________________
//

VarLabel*
VarLabel::create(const std::string& name, const TypeDescription* td,
                 const IntVector& boundaryLayer, /* = IntVector(0,0,0) */
                 VarType vartype                 /* = Normal */
)
{
  VarLabel* label;

  g_label_mutex.lock();
  {
    auto iter = g_all_labels.find(name);
    if (iter != g_all_labels.end()) {
      // two labels with the same name -- make sure they are the same type
      VarLabel* dup = iter->second;
      if (boundaryLayer != dup->d_boundary_layer) {
        std::ostringstream out;
        out << "Multiple VarLabels for " << dup->getName()
            << " defined with different # of boundary layers";
        SCI_THROW(InternalError(out.str(), __FILE__, __LINE__));
      }

      if (td != dup->d_td || vartype != dup->d_var_type) {
        std::ostringstream out;
        out << "VarLabel with same name exists, '" << name
            << "', but with different type";
        SCI_THROW(InternalError(out.str(), __FILE__, __LINE__));
      }

      label = dup;
    } else {
      label = scinew VarLabel(name, td, boundaryLayer, vartype);
      g_all_labels[name] = label;
      //DEBUGOUT(g_varlabel_dbg,
      //         "Created VarLabel: " << label->d_name << " [address = " << label);
      dbg << "Created VarLabel: " << label->d_name << " [address = " << label
          << "\n";
    }
    label->addReference();
  }
  g_label_mutex.unlock();

  return label;
}

bool
VarLabel::destroy(const VarLabel* label)
{
  if (label == nullptr) {
    return false;
  }

  if (label->removeReference()) {
    g_label_mutex.lock();
    {
      auto iter = g_all_labels.find(label->d_name);
      if (iter != g_all_labels.end() && iter->second == label) {
        g_all_labels.erase(iter);
      }
      //DEBUGOUT(g_varlabel_dbg,
      //           "Deleted VarLabel: " << label->d_name);
      dbg << "Deleted VarLabel: " << label->d_name << std::endl;
    }
    g_label_mutex.unlock();

    delete label;

    return true;
  }

  return false;
}

VarLabel::VarLabel(const std::string& name, const Uintah::TypeDescription* td,
                   const IntVector& boundaryLayer, VarType vartype)
  : d_name(name)
  , d_td(td)
  , d_boundary_layer(boundaryLayer)
  , d_var_type(vartype)
{
}

void
VarLabel::printAll()
{
  for (auto label : g_all_labels) {
    std::cout << label.second->d_name << std::endl;
  }
  //auto iter = g_all_labels.begin();

  //for (; iter != g_all_labels.end(); iter++) {
  //  std::cout << (*iter).second->d_name << std::endl;
  //}
}

VarLabel*
VarLabel::find(const std::string& name)
{
  auto found = g_all_labels.find(name);

  if (found == g_all_labels.end()) {
    return nullptr;
  } else {
    return found->second;
  }
}

VarLabel*
VarLabel::find(const std::string& name,
               const std::string& message)
{
  auto label = VarLabel::find(name);
  if (label == nullptr) {
    std::ostringstream warn;
    warn << message;
    warn << "**ERROR** Could not find the VarLabel (" << name << " ).";
    throw InternalError(warn.str(), __FILE__, __LINE__);
    return nullptr;
  }

  return label;
}

VarLabel*
VarLabel::particlePositionLabel()
{
  return find(s_particle_position_name);
}

std::string
VarLabel::getFullName(int matlIndex, const Patch* patch) const
{
  std::ostringstream out;
  out << d_name << "(matl=" << matlIndex;

  if (patch) {
    out << ", patch=" << patch->getID();
  } else {
    out << ", no patch";
  }
  out << ")";

  return out.str();
}

void
VarLabel::isReductionTask(bool input)
{
  if (!d_td->isReductionVariable()) {
    std::ostringstream out;
    out << "Only reduction variables may allow multiple computes.\n'" 
        << d_name << "' is not a reduction variable.";
    SCI_THROW(InternalError(out.str(), __FILE__, __LINE__));
  }

  d_is_reduction_task = input;
}

std::ostream&
operator<<(std::ostream& out, const Uintah::VarLabel& vl)
{
  out << vl.getName();
  return out;
}

} // namespace Uintah
