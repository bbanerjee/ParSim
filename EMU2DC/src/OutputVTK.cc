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
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkSmartPointer.h>

#include <unistd.h>
#include <fstream>
#include <regex>

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
  std::string name_without_ext = outputFile();
  int lastIndex = name_without_ext.find_last_of(".");
  if (lastIndex != std::string::npos) {
    name_without_ext = name_without_ext.substr(0, lastIndex); 
  }
  
  std::ostringstream of_directory;
  of_directory << "./" << name_without_ext <<"/";
  std::ostringstream of_name;
  of_name << of_directory.str() << name_without_ext << std::setfill('0') << std::setw(5) << outputFileCount() << ".vtu"; 

  // Create a writer
  vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer = 
     vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
  writer->SetFileName((of_name.str()).c_str());
  //writer->SetFileName(name_without_ext.c_str());
  writer->SetTimeStep(outputFileCount());

  // Create a pointer to a VTK MultiBlock data set
  vtkSmartPointer<vtkMultiBlockDataSet> data_set = vtkSmartPointer<vtkMultiBlockDataSet>::New();

  // Loop through bodies
  int body_count = 1;
  for (auto body_iter = bodyList.begin(); body_iter != bodyList.end(); ++body_iter) {

    // Get the node list for the body
    const NodePArray& node_list = (*body_iter)->nodes();

    // Create pointer to VTK UnstructuredGrid data set
    vtkSmartPointer<vtkUnstructuredGrid> point_data = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Add the time
    addTimeToVTKDataSet(time.currentTime(), point_data);

    // Create the unstructured data set for the body
    createVTKUnstructuredDataSet(node_list, point_data);

    // Add the unstructured data to the multi block
    data_set->SetBlock(body_count, point_data);

    // Increment body count
    ++body_count;
  }

  // Write the data
  writer->SetInput(data_set);
  writer->Write();

  // Increment the output file count
  incrementOutputFileCount();
}

void
OutputVTK::createVTKUnstructuredDataSet(const NodePArray& nodeList, 
                                        vtkSmartPointer<vtkUnstructuredGrid>& dataSet)
{

  // Find number of points
  int count = 1;
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
  damage->SetNumberOfComponents(1);
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
    //std::cout << "Damage array = " << damage << std::endl;
    //std::cout << "size = " << nodeList.size() << "count = " << count << " id = " << id << " index = " << cur_node->damageIndex() << std::endl;
    damage->InsertValue(id, cur_node->damageIndex());
    disp->InsertTuple(id, displacement);
    vel->InsertTuple(id, velocity);
  }

  // Add points to data set
  dataSet->SetPoints(pts);
  dataSet->GetPointData()->AddArray(damage);
  dataSet->GetPointData()->AddArray(disp);
  dataSet->GetPointData()->AddArray(vel);

  // Check point data
  vtkPointData *pd = dataSet->GetPointData();
  if (pd) {
    std::cout << " contains point data with " << pd->GetNumberOfArrays() << " arrays." << std::endl;
    for (int i = 0; i < pd->GetNumberOfArrays(); i++) {
      std::cout << "\tArray " << i << " is named "
                << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL")
                << std::endl;
    }
  }
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
