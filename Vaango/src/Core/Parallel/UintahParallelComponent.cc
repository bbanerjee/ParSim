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

#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/UintahParallelComponent.h>

#include <algorithm>

namespace Uintah {

UintahParallelComponent::UintahParallelComponent(const ProcessorGroup* myworld)
  : d_myworld(myworld)
{
}

UintahParallelComponent::~UintahParallelComponent() noexcept(false) {}

void
UintahParallelComponent::attachPort(const std::string& name,
                                    UintahParallelPort* port)
{
  auto iter = portmap.find(name);
  if (iter == portmap.end()) {
    portmap[name] = std::make_unique<PortRecord>(port);
  } else {
    iter->second->connections.push_back(port);
  }
}

UintahParallelComponent::PortRecord::PortRecord(UintahParallelPort* port)
{
  connections.push_back(port);
}

UintahParallelPort*
UintahParallelComponent::getPort(const std::string& name)
{
  auto iter = portmap.find(name);
  if (iter == portmap.end()) {
    return nullptr;
  } else if (iter->second->connections.size() > 1) {
    return iter->second->connections.back();
  } else {
    return iter->second->connections[0];
  }
}

UintahParallelPort*
UintahParallelComponent::getPort(const std::string& name, unsigned int i)
{
  auto iter = portmap.find(name);
  if (iter == portmap.end()) {
    return 0;
  } else if (iter->second->connections.size() > 1) {
    return iter->second->connections[i];
  } else {
    return iter->second->connections[0];
  }
}

void
UintahParallelComponent::releasePort(const std::string&)
{
}

unsigned int
UintahParallelComponent::numConnections(const std::string& name)
{
  auto iter = portmap.find(name);
  if (iter == portmap.end()) {
    return 0;
  } else {
    return iter->second->connections.size();
  }
}

} // end namespace Uintah