/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#ifndef __MATITI_OUTPUT_VTK_H__
#define __MATITI_OUTPUT_VTK_H__

#include <InputOutput/Output.h>
#include <Containers/NodePArray.h>

#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>

#include <iostream>
#include <sstream>

namespace Matiti {

  class OutputVTK : public Output
  {
  public:
    OutputVTK();
    OutputVTK(const std::string& fileName,
              int iterInterval);
    OutputVTK(const Uintah::ProblemSpecP& ps);
    virtual ~OutputVTK();


    void write(const Time& time, const Domain& domain, const BodySPArray& bodyList);
    void write(const Time& time, const Domain& domain, const RigidBodySPArray& bodyList);

    void writeMB(const Time& time, const Domain& domain, const BodySPArray& bodyList);

  private:

    void getFileNames(std::ostringstream& domainFileName,
                      std::ostringstream& nodeFileName);

    void writeDomain(const Time& time, const Domain& domain, 
                     std::ostringstream& fileName);

    void writeNodes(const Time& time, const BodySPArray& bodyList,
                    std::ostringstream& fileName);
 
    void createVTKUnstructuredGrid(const NodePArray& nodeList, 
                                   vtkSmartPointer<vtkPoints>& pts,
                                   vtkSmartPointer<vtkUnstructuredGrid>& dataSet);

    void createVTKUnstructuredDataSet(const Domain& domain,
                                      const NodePArray& nodeList, 
                                      vtkSmartPointer<vtkUnstructuredGrid>& dataSet);

    void addTimeToVTKDataSet(double time, 
                             vtkSmartPointer<vtkUnstructuredGrid>& dataSet);

    void addDomainToVTKDataSet(const Domain& domain, 
                               vtkSmartPointer<vtkUnstructuredGrid>& dataSet);

    void addDomainToVTKUnstructuredGrid(const Domain& domain, 
                                        vtkSmartPointer<vtkPoints>& pts,
                                        vtkSmartPointer<vtkUnstructuredGrid>& dataSet);

  private:

    std::ostringstream d_output_dir;

  }; // end class

} // end namespace

#endif

