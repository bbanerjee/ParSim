#include <Boundary/PlaneBoundary.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(PlaneBoundaryTest, readCSV) {

  Boundary::BoundaryType boundary_type = Boundary::BoundaryType::PLANE;
  int extraNum = 0;
  int id = 1;
  std::vector<double> values = {{-1, 0, 0, 0, 0, 25}};

  // Create csv data 
  std::stringstream csv;
  csv << extraNum << "\n";
  csv << id;
  for (int ii = 0; ii < 6; ii++) {
    csv << values[ii] << " ";
  }
  csv << "\n";

  // Read csv data
  std::istream& stream(csv);
  PlaneBoundary plane_csv(boundary_type, stream);
  //std::cout << "PlaneBoundary: CSV:" << std::endl;
  //std::cout << plane_csv;

  // Add traction BC to csv data
  std::stringstream csv_traction_bc;
  csv_traction_bc << extraNum << "\n";
  csv_traction_bc << id;
  for (int ii = 0; ii < 6; ii++) {
    csv_traction_bc << values[ii] << " ";
  }
  csv_traction_bc << "\n";
  int bcType = static_cast<int>(PlaneBoundary::BCType::TRACTION);
  csv_traction_bc << bcType << "\n";

  int num_points = 3;
  int num_components = 1;
  std::vector<double> times(num_points);
  std::vector<double> bc_values(num_points);
  for (int ii = 0; ii < num_points; ii++) {
    times[ii] = ii+1;
    bc_values[ii] = 5*(ii+1);
  }

  csv_traction_bc << num_points << " " << num_components << "\n";
  for (int ii = 0; ii < num_points; ii++) {
    csv_traction_bc << times[ii] << " ";
  }
  csv_traction_bc << "\n";
  for (int ii = 0; ii < num_points; ii++) {
    csv_traction_bc << bc_values[ii] << " ";
  }
  csv_traction_bc << "\n";
  //std::cout << csv_traction_bc.str();

  // Read csv data
  std::istream& stream_traction_bc(csv_traction_bc);
  PlaneBoundary plane_csv_traction_bc(boundary_type, stream_traction_bc);
  //std::cout << "PlaneBoundary: CSV with traction BC:" << std::endl;
  //std::cout << plane_csv_traction_bc;

  // Add displacement BC to csv data
  std::stringstream csv_disp_bc;
  csv_disp_bc << extraNum << "\n";
  csv_disp_bc << id;
  for (int ii = 0; ii < 6; ii++) {
    csv_disp_bc << values[ii] << " ";
  }
  csv_disp_bc << "\n";
  bcType = static_cast<int>(PlaneBoundary::BCType::DISPLACEMENT);
  csv_disp_bc << bcType << "\n";

  csv_disp_bc << num_points << " " << num_components << "\n";
  for (int ii = 0; ii < num_points; ii++) {
    csv_disp_bc << times[ii] << " ";
  }
  csv_disp_bc << "\n";
  for (int ii = 0; ii < num_points; ii++) {
    csv_disp_bc << bc_values[ii] << " ";
  }
  csv_disp_bc << "\n";
  //std::cout << csv_disp_bc.str();

  // Read csv data
  std::istream& stream_disp_bc(csv_disp_bc);
  PlaneBoundary plane_csv_disp_bc(boundary_type, stream_disp_bc);
  //std::cout << "PlaneBoundary: CSV with disp BC:" << std::endl;
  //std::cout << plane_csv_disp_bc;
}

TEST(PlaneBoundaryTest, readXML) {

  Boundary::BoundaryType boundary_type = Boundary::BoundaryType::PLANE;
  Boundary::BoundaryID boundary_id = static_cast<Boundary::BoundaryID>(0);

  zen::XmlDoc doc;
  zen::XmlOut xml(doc);
  xml["direction"]("[-1, 0, 0]");
  xml["position"]("[0, 0, 25]");
  xml["initial_velocity"]("[10, 0, 0]");
  zen::XmlIn xml_in(doc);

  PlaneBoundary plane_xml(boundary_type, boundary_id, xml_in);
  //std::cout << "PlaneBoundary: XML with no BC:" << std::endl;
  //std::cout << plane_xml;

  // With traction bc
  int num_points = 3;
  std::vector<double> times(num_points);
  std::vector<double> bc_values(num_points);
  for (int ii = 0; ii < num_points; ii++) {
    times[ii] = ii+1;
    bc_values[ii] = 5*(ii+1);
  }

  zen::XmlDoc doc_traction;
  zen::XmlOut xml_traction(doc_traction);
  xml_traction["direction"]("[-1, 0, 0]");
  xml_traction["position"]("[0, 0, 25]");
  xml_traction["initial_velocity"]("[10, 0, 0]");
  std::string timeStr = "";
  std::string valueStr = "";
  for (int ii = 0; ii < num_points; ii++) {
    timeStr += std::to_string(times[ii]) + " ";
    valueStr += std::to_string(bc_values[ii]) + " ";
  }
  xml_traction["tractionBC"]["time"](timeStr);
  xml_traction["tractionBC"]["value"](valueStr);
  zen::XmlIn xml_in_traction(doc_traction);

  PlaneBoundary plane_xml_traction(boundary_type, boundary_id, xml_in_traction);
  //std::cout << "PlaneBoundary: XML with traction BC:" << std::endl;
  //std::cout << plane_xml_traction;

  // With displacement bc
  zen::XmlDoc doc_disp;
  zen::XmlOut xml_disp(doc_disp);
  xml_disp["direction"]("[-1, 0, 0]");
  xml_disp["position"]("[0, 0, 25]");
  xml_disp["initial_velocity"]("[10, 0, 0]");
  xml_disp["displacementBC"]["time"](timeStr);
  xml_disp["displacementBC"]["value"](valueStr);
  zen::XmlIn xml_in_disp(doc_disp);
  //std::cout << zen::serialize(doc_disp);

  PlaneBoundary plane_xml_disp(boundary_type, boundary_id, xml_in_disp);
  //std::cout << "PlaneBoundary: XML with displacement BC:" << std::endl;
  //std::cout << plane_xml_disp;
}

TEST(PlaneBoundaryTest, readJSON) {

  Boundary::BoundaryType boundary_type = Boundary::BoundaryType::PLANE;
  Boundary::BoundaryID boundary_id = static_cast<Boundary::BoundaryID>(0);

  nlohmann::json json_data;
  json_data["direction"] = "[-1, 0, 0]";
  json_data["position"] = "[0, 0, 25]";
  json_data["initial_velocity"] = "[10, 0, 0]";

  PlaneBoundary plane_json(boundary_type, boundary_id, json_data);
  //std::cout << plane_json;

  // With traction bc
  int num_points = 3;
  std::vector<double> times(num_points);
  std::vector<double> bc_values(num_points);
  for (int ii = 0; ii < num_points; ii++) {
    times[ii] = ii+1;
    bc_values[ii] = 5*(ii+1);
  }
  std::string timeStr = "";
  std::string valueStr = "";
  for (int ii = 0; ii < num_points; ii++) {
    timeStr += std::to_string(times[ii]) + " ";
    valueStr += std::to_string(bc_values[ii]) + " ";
  }
  nlohmann::json json_data_traction = json_data;
  json_data_traction["tractionBC"]["time"] = timeStr;
  json_data_traction["tractionBC"]["value"] = valueStr;
  //std::cout << json_data_traction;

  PlaneBoundary plane_json_traction(boundary_type, boundary_id, json_data_traction);
  //std::cout << "PlaneBoundary: JSON with traction BC:" << std::endl;
  //std::cout << plane_json_traction;

  // With displacement bc
  nlohmann::json json_data_disp = json_data;
  json_data_disp["displacementBC"]["time"] = timeStr;
  json_data_disp["displacementBC"]["value"] = valueStr;
  //std::cout << json_data_disp;

  PlaneBoundary plane_json_disp(boundary_type, boundary_id, json_data_disp);
  //std::cout << "PlaneBoundary: JSON with disp BC:" << std::endl;
  //std::cout << plane_json_disp;
}
