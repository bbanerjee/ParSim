#include <OutputVTK.h>
#include <Exception.h>
#include <NodePArray.h>
#include <Body.h>
#include <Node.h>

#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkSmartPointer.h>

#include <unistd.h>
#include <fstream>

using namespace Emu2DC;

OutputVTK::OutputVTK()
  : Output()
{
}

OutputVTK::OutputVTK(const Uintah::ProblemSpecP& ps)
  : Output(ps)
{
}

OutputVTK::~OutputVTK()
{
}

void
OutputVTK::write(const Time& time, const BodySPArray& bodyList) 
{
  // Write the output to individual files
  std::ostringstream of_name;
  of_name << outputFile() << std::setw(5) << outputFileCount() << ".vtu"; 

  // Create a pointer to a VTK MultiBlock data set
  vtkSmartPointer<vtkMultiBlockDataSet> data_set = vtkSmartPointer<vtkMultiBlockDataSet>::New();

  // Loop through bodies
  int body_count = 0;
  for (auto body_iter = bodyList.begin(); body_iter != bodyList.end(); ++body_iter) {

    // Get the node list for the body
    const NodePArray& node_list = (*body_iter)->nodes();

    // Create pointer to VTK UnstructuredGrid data set
    //vtkSmartPointer<vtkUnstructuredGrid> point_data = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkUnstructuredGrid> point_data = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Add the time
    addTimeToVTKDataSet(time.currentTime(), point_data);

    // Create the unstructured data set for the body
    createVTKUnstructuredDataSet(node_list, point_data);

    // Add the unstructured data to the multi block
    data_set->SetBlock(body_count, point_data);

    // Delete
    point_data->Delete();
    ++body_count;
  }

  // Write the data
  vtkSmartPointer<vtkGenericDataObjectWriter> writer = 
     vtkSmartPointer<vtkGenericDataObjectWriter>::New();
  writer->SetInput(data_set);
  writer->SetFileName((of_name.str()).c_str());
  writer->Write();

  // Increment the output file count
  incrementOutputFileCount();
}

void
OutputVTK::createVTKUnstructuredDataSet(const NodePArray& nodeList, 
                                        vtkSmartPointer<vtkUnstructuredGrid>& dataSet)
{

  // Find number of points
  int count = 0;
  for (auto node_iter = nodeList.begin(); node_iter != nodeList.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    if (cur_node->omit()) continue;  // skip this node
    ++count;
  }

  // Set up pointer to point data
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New(); 
  pts->SetNumberOfPoints(count);
  
  // Set up pointer for damage data
  vtkSmartPointer<vtkDoubleArray> damage = vtkSmartPointer<vtkDoubleArray>::New();
  damage->SetNumberOfTuples(count);
  damage->SetName("Damage");

  // Set up pointer for displacement and velocity data
  vtkSmartPointer<vtkDoubleArray> disp = vtkSmartPointer<vtkDoubleArray>::New();
  disp->SetNumberOfComponents(3);
  disp->SetNumberOfTuples(count);
  disp->SetName("Displacement");

  vtkSmartPointer<vtkDoubleArray> vel = vtkSmartPointer<vtkDoubleArray>::New();
  vel->SetNumberOfComponents(3);
  vel->SetNumberOfTuples(count);
  vel->SetName("Velocity");

  // Loop through nodes
  int id = 0;
  for (auto node_iter = nodeList.begin(); node_iter != nodeList.end(); ++node_iter) {
    NodeP cur_node = *node_iter;
    if (cur_node->omit()) continue;  // skip this node
    ++id;
    double displacement[3], position[3], velocity[3];
    for (int ii = 0; ii < 3; ++ii) {
      displacement[ii] = cur_node->displacement()[ii];
      position[ii] = cur_node->position()[ii] + displacement[ii];
      velocity[ii] = cur_node->velocity()[ii];
    }
    pts->SetPoint(id, position);
    damage->InsertValue(id, cur_node->damageIndex());
    disp->InsertTuple(id, displacement);
    vel->InsertTuple(id, velocity);
  }

  // Add points to data set
  dataSet->SetPoints(pts);
  dataSet->GetPointData()->SetScalars(damage);
  dataSet->GetPointData()->SetVectors(disp);
  dataSet->GetPointData()->SetVectors(vel);
}

void 
OutputVTK::addTimeToVTKDataSet(double time, vtkSmartPointer<vtkUnstructuredGrid>& dataSet)
{
  vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetName("TIME");
  array->SetNumberOfTuples(1);
  array->SetTuple1(0, time);
  dataSet->GetFieldData()->AddArray(array);
}
