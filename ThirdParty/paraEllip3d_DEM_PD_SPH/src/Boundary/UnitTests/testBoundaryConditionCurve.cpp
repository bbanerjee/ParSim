#include <Boundary/BoundaryConditionCurve.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(BoundaryConditionCurveTest, test1D) {

  int num_points = 3;
  int num_components = 1;

  std::vector<double> times(num_points);
  std::vector<double> values(num_points);
  for (int ii = 0; ii < num_points; ii++) {
    times[ii] = ii+1;
    values[ii] = 5*(ii+1);
  }

  /* Read csv data */
  std::stringstream csv;
  csv << num_points << " " << num_components << "\n";
  for (int ii = 0; ii < num_points; ii++) {
    csv << times[ii] << " ";
  }
  csv << "\n";

  for (int ii = 0; ii < num_points; ii++) {
    csv << values[ii] << " ";
  }
  csv << "\n";

  std::istream& stream(csv);
  BoundaryConditionCurve<double> bc_csv(stream);
  //std::cout << "BC: 1D : CSV:" << std::endl;
  //std::cout << bc_csv;

  /* Read xml data */
  zen::XmlDoc doc;
  zen::XmlOut xml(doc);
  std::string timeStr = "";
  std::string valueStr = "";
  for (int ii = 0; ii < num_points; ii++) {
    timeStr += std::to_string(times[ii]) + " ";
    valueStr += std::to_string(values[ii]) + " ";
  }
  xml["time"](timeStr);
  xml["value"](valueStr);
  zen::XmlIn inxml(doc);
  BoundaryConditionCurve<double> bc_xml(inxml);
  //std::cout << "BC: 1D : XML:" << std::endl;
  //std::cout << bc_xml;

  /* Read json data */
  nlohmann::json json_data;
  json_data["time"] = timeStr;
  json_data["value"] = valueStr;
  BoundaryConditionCurve<double> bc_json(json_data);
  //std::cout << "BC: 1D : JSON:" << std::endl;
  //std::cout << bc_json;

  /* interpolate */
  double val = bc_xml.getBCValue(0);
  EXPECT_EQ(val, 5);
  //std::cout << "val = " << val << "\n";

  val = bc_xml.getBCValue(5);
  EXPECT_EQ(val, 15);
  //std::cout << "val = " << val << "\n";

  val = bc_xml.getBCValue(2.999999999999);
  EXPECT_NEAR(val, 15, 1.0e-8);
  //std::cout << "val = " << val << "\n";

  val = bc_xml.getBCValue(1.37);
  double ss = (1.37 - 1)/(2 - 1);
  double vv = (1 - ss)*5 + ss*10;
  EXPECT_EQ(val, vv);
  //std::cout << "val = " << val << "\n";
}

TEST(BoundaryConditionCurveTest, test3D) {

  int num_points = 3;
  int num_components = 3;

  std::vector<double> times(num_points);
  std::vector<Vec> values(num_points);
  for (int ii = 0; ii < num_points; ii++) {
    times[ii] = ii+1;
    values[ii] = Vec(5*(ii+1), 6*(ii+1), 7*(ii+1));
  }

  /* Read csv data */
  std::stringstream csv;
  csv << num_points << " " << num_components << "\n";
  for (int ii = 0; ii < num_points; ii++) {
    csv << times[ii] << " ";
  }
  csv << "\n";

  for (int ii = 0; ii < num_points; ii++) {
    csv << values[ii] << " ";
  }
  csv << "\n";

  std::istream& stream(csv);
  BoundaryConditionCurve<Vec> bc_csv(stream);
  //std::cout << "BC: 3D : CSV:" << std::endl;
  //std::cout << bc_csv;

  /* Read xml data */
  zen::XmlDoc doc;
  zen::XmlOut xml(doc);
  std::string timeStr = "";
  std::string valueStr = "";
  for (int ii = 0; ii < num_points; ii++) {
    timeStr += std::to_string(times[ii]) + " ";
    valueStr += std::to_string(values[ii].x()) + " " +
                std::to_string(values[ii].y()) + " " +
                std::to_string(values[ii].z()) + " ";
  }
  xml["time"](timeStr);
  xml["value"](valueStr);
  zen::XmlIn inxml(doc);
  BoundaryConditionCurve<Vec> bc_xml(inxml);
  //std::cout << "BC: 3D : XML:" << std::endl;
  //std::cout << bc_xml;

  /* Read json data */
  nlohmann::json json_data;
  json_data["time"] = timeStr;
  json_data["value"] = valueStr;
  BoundaryConditionCurve<Vec> bc_json(json_data);
  //std::cout << "BC: 3D : JSON:" << std::endl;
  //std::cout << bc_json;

  /* interpolate */
  Vec val = bc_xml.getBCValue(0);
  EXPECT_EQ(val.x(), 5);
  EXPECT_EQ(val.y(), 6);
  EXPECT_EQ(val.z(), 7);
  //std::cout << "val = " << val << "\n";

  val = bc_xml.getBCValue(5);
  //std::cout << "val = " << val << "\n";
  EXPECT_EQ(val.x(), 15);
  EXPECT_EQ(val.y(), 18);
  EXPECT_EQ(val.z(), 21);

  val = bc_xml.getBCValue(2.999999999999);
  //std::cout << "val = " << val << "\n";
  EXPECT_NEAR(val.x(), 15, 1.0e-8);
  EXPECT_NEAR(val.y(), 18, 1.0e-8);
  EXPECT_NEAR(val.z(), 21, 1.0e-8);

  val = bc_xml.getBCValue(1.37);
  double ss = (1.37 - 1)/(2 - 1);
  Vec vv = (1 - ss)*Vec(5,6,7) + ss*Vec(10,12,14);
  //std::cout << "val = " << val << " vv = " << vv << "\n";
  EXPECT_EQ(val.x(), vv.x());
  EXPECT_EQ(val.y(), vv.y());
  EXPECT_EQ(val.z(), vv.z());
}