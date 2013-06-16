#ifndef __EMU2DC_OUTPUT_VTK_H__
#define __EMU2DC_OUTPUT_VTK_H__

#include <Output.h>
#include <NodePArray.h>

#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>

#include <iostream>

namespace Emu2DC {

  class OutputVTK : public Output
  {
  public:
    OutputVTK();
    OutputVTK(const Uintah::ProblemSpecP& ps);
    virtual ~OutputVTK();

    void write(const Time& time, const Domain& domain, const BodySPArray& bodyList);

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

  }; // end class

} // end namespace

#endif

