/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#ifndef __CCA_COMPONENTS_LOADBALANCERS_RoundRobinLoadBalancer_H__
#define __CCA_COMPONENTS_LOADBALANCERS_RoundRobinLoadBalancer_H__

#include <CCA/Components/LoadBalancers/LoadBalancerCommon.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <set>

namespace Uintah {

class RoundRobinLoadBalancer : public LoadBalancerCommon
{
public:
  RoundRobinLoadBalancer(const ProcessorGroup* myworld);
  ~RoundRobinLoadBalancer() = default;

  RoundRobinLoadBalancer(const RoundRobinLoadBalancer&) = delete;
  RoundRobinLoadBalancer(RoundRobinLoadBalancer&&)      = delete;
  RoundRobinLoadBalancer&
  operator=(const RoundRobinLoadBalancer&) = delete;
  RoundRobinLoadBalancer&
  operator=(RoundRobinLoadBalancer&&) = delete;

  virtual int
  getPatchwiseProcessorAssignment(const Patch* patch);

  virtual bool
  needRecompile(const GridP&)
  {
    return false;
  };

private:
};
} // End namespace Uintah

#endif //__CCA_COMPONENTS_LOADBALANCERS_RoundRobinLoadBalancer_H__
