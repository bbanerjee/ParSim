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

#include <Core/Disclosure/TypeDescription.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/Assert.h>

//#include <Core/Parallel/CrowdMonitor.h>
// 

#include <Core/Parallel/CrowdMonitor.h>
#include <Core/Parallel/MasterLock.h>

#include <iostream>
#include <sstream>
#include <map>
#include <vector>

namespace {

// static Mutex tdLock("TypeDescription::getMPIType lock");
// static CrowdMonitor tpLock("TypeDescription type lock");

using register_monitor =
  Uintah::CrowdMonitor<Uintah::TypeDescription::register_tag>;
using lookup_monitor =
  Uintah::CrowdMonitor<Uintah::TypeDescription::lookup_tag>;

std::map<std::string, const Uintah::TypeDescription*>* g_types = nullptr;
std::vector<const Uintah::TypeDescription*>* g_typelist = nullptr;
bool g_killed = false;

Uintah::MasterLock get_mpi_type_lock{};
} // namespace

namespace Uintah {

TypeDescription::TypeDescription(Type type, const std::string& name,
                                 bool isFlat, MPI_Datatype (*mpitypemaker)())
  : d_type{ type }
  , d_name{ name }
  , d_isFlat{ isFlat }
  , d_mpitype{ MPI_Datatype(-1) }
  , d_mpitypemaker{ mpitypemaker }
{
  register_type();
}

TypeDescription::TypeDescription(Type type, const std::string& name,
                                 bool isFlat, MPI_Datatype mpitype)
  : d_type{ type }
  , d_name{ name }
  , d_isFlat{ isFlat }
  , d_mpitype{ mpitype }
{
  register_type();
}

TypeDescription::TypeDescription(Type type, const std::string& name,
                                 Variable* (*maker)(),
                                 const TypeDescription* subtype)
  : d_type{ type }
  , d_subtype(subtype)
  , d_name{ name }
  , d_mpitype{ MPI_Datatype(-2) }
  , d_maker{ maker }
{
  register_type();
}

void
TypeDescription::register_type()
{
  {
    register_monitor register_write_lock{
      Uintah::CrowdMonitor<TypeDescription::register_tag>::WRITER
    };

    if (!g_types) {
      ASSERT(!g_killed);
      ASSERT(!g_typelist)
      g_types = scinew std::map<std::string, const TypeDescription*>;
      g_typelist = scinew std::vector<const TypeDescription*>;
    }

    auto iter = g_types->find(getName());
    if (iter == g_types->end()) {
      (*g_types)[getName()] = this;
    }
    g_typelist->push_back(this);
  }
}

void
TypeDescription::deleteAll()
{
  if (!g_types) {
    ASSERT(!g_killed);
    ASSERT(!g_typelist);
    return;
  }
  g_killed = true;
  for (auto iter = g_typelist->begin(); iter != g_typelist->end(); iter++) {
    delete *iter;
  }
  delete g_types;
  delete g_typelist;

  g_types = nullptr;
  g_typelist = nullptr;
}

std::string
TypeDescription::getName() const
{
  if (d_subtype) {
    return d_name + "<" + d_subtype->getName() + ">";
  } else {
    return d_name;
  }
}

std::string
TypeDescription::getFileName() const
{
  if (d_subtype) {
    return d_name + d_subtype->getFileName();
  } else {
    return d_name;
  }
}

const TypeDescription*
TypeDescription::lookupType(const std::string& t)
{
  {
    lookup_monitor lookup_read_lock{
      Uintah::CrowdMonitor<TypeDescription::lookup_tag>::READER
    };

    if (!g_types) {
      return 0;
    }

    auto iter = g_types->find(t);
    if (iter == g_types->end()) {
      return 0;
    }

    return iter->second;
  }
}

MPI_Datatype
TypeDescription::getMPIType() const
{
  if (d_mpitype == MPI_Datatype(-1)) {
    {
      std::lock_guard<Uintah::MasterLock> guard(get_mpi_type_lock);
      if (d_mpitype == MPI_Datatype(-1)) {
        if (d_mpitypemaker) {
          d_mpitype = (*d_mpitypemaker)();
        } else {
          std::ostringstream out;
          out << "MPI Datatype requested, but do not know how to make it";
          throw InternalError(out.str(), __FILE__, __LINE__);
        }
      }
    }
  }
  ASSERT(d_mpitype != MPI_Datatype(-2));
  return d_mpitype;
}

Variable*
TypeDescription::createInstance() const
{
  if (!d_maker) {
    std::ostringstream out;
    out << "Do not know how to create instance for type: " << getName(),
      throw InternalError(out.str(), __FILE__, __LINE__);
  }
  return (*d_maker)();
}

std::string
TypeDescription::toString(const Type& type)
{
  switch (type) {
    case Type::CCVariable:
      return "CCVariable";
    case Type::NCVariable:
      return "NCVariable";
    case Type::SFCXVariable:
      return "SFCXVariable";
    case Type::SFCYVariable:
      return "SFCYVariable";
    case Type::SFCZVariable:
      return "SFCZVariable";
    case Type::ParticleVariable:
      return "ParticleVariable";
    case Type::PerPatch:
      return "PerPatch";
    case Type::Point:
      return "Point";
    case Type::Vector:
      return "Vector";
    case Type::Matrix3:
      return "Matrix3";
    case Type::ReductionVariable:
      return "ReductionVariable";
    case Type::SoleVariable:
      return "SoleVariable";
    case Type::double_type:
      return "double_type";
    case Type::float_type:
      return "float_type";
    case Type::bool_type:
      return "bool_type";
    case Type::int_type:
      return "int_type";
    case Type::short_int_type:
      return "short_int_type";
    case Type::long_type:
      return "long_type";
    case Type::long64_type:
      return "long64_type";
    case Type::Short27:
      return "Short27";
    case Type::Int130:
      return "Int130";
    case Type::Stencil4:
      return "Stencil4";
    case Type::Stencil7:
      return "Stencil7";
    case Type::NeighborList:              // for Peridynamics
      return "NeighborList";
    case Type::NeighborConnectivity:      // for Peridynamics
      return "NeighborConnectivity";
    case Type::NeighborBondEnergy:        // for Peridynamics
      return "NeighborBondEnergy";
    case Type::NeighborBondInternalForce: // for Peridynamics
      return "NeighborBondInternalForce";
    case Type::MetalIntvar:
      return "MetalIntVar";
    case Type::ArenaIntvar:
      return "ArenaIntVar";
    case Type::BorjaIntvar:
      return "BorjaIntVar";
    case Type::SoilBrannonIntvar:
      return "SoilBrannonIntVar";
    case Type::TabularCapIntvar:
      return "TabularCapIntVar";
    case Type::IntVector:
      return "IntVector";
    case Type::Unknown:
      return "Unknown";
    case Type::Other:
      return "Other";
    default:
      std::stringstream msg;
      msg << "Invalid type: " << static_cast<std::underlying_type<Type>::type>(type);
      throw InternalError(msg.str(), __FILE__, __LINE__);
  }
}
} // end namespace Uintah
