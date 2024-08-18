#ifndef __CORE_PARALLEL_UNITAHPARALLELCOMPONENT_H__
#define __CORE_PARALLEL_UNITAHPARALLELCOMPONENT_H__

/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace Uintah {

class UintahParallelPort;
class ProcessorGroup;

class UintahParallelComponent
{
  struct PortRecord
  {
    PortRecord(UintahParallelPort* conn);
    std::vector<UintahParallelPort*> connections;
  };
  std::map<std::string, std::unique_ptr<PortRecord>> portmap;

public:
  UintahParallelComponent(const ProcessorGroup* myworld);
  virtual ~UintahParallelComponent() noexcept(false);

  void
  attachPort(const std::string& name, UintahParallelPort* port);

  UintahParallelPort*
  getPort(const std::string& name);

  UintahParallelPort*
  getPort(const std::string& name, unsigned int i);

  void
  releasePort(const std::string& name);

  unsigned int
  numConnections(const std::string& name);

  virtual void
  setComponents(UintahParallelComponent* comp) = 0;

  virtual void
  getComponents() = 0;

  virtual void
  releaseComponents() = 0;

protected:
  const ProcessorGroup* d_myworld;
};
} // End namespace Uintah

#endif //__CORE_PARALLEL_UNITAHPARALLELCOMPONENT_H__
