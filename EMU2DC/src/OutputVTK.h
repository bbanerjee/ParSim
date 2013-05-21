#ifndef __EMU2DC_OUTPUT_VTK_H__
#define __EMU2DC_OUTPUT_VTK_H__

#include <Output.h>
#include <NodePArray.h>

#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

namespace Emu2DC {

  class OutputVTK : public Output
  {
  public:
    OutputVTK();
    OutputVTK(const Uintah::ProblemSpecP& ps);
    virtual ~OutputVTK();

    void write(const Time& time, const BodySPArray& bodyList);

  private:
    void createVTKUnstructuredDataSet(const NodePArray& nodeList, 
                                      vtkSmartPointer<vtkUnstructuredGrid>& dataSet);

    void addTimeToVTKDataSet(double time, 
                             vtkSmartPointer<vtkUnstructuredGrid>& dataSet);

    

  }; // end class

} // end namespace

#endif

