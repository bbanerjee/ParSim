#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_UTIL_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_UTIL_H

#include <string>
#include <vector>

namespace Vaango {

namespace Util {

// Taken from:
// https://stackoverflow.com/questions/236129/most-elegant-way-to-split-a-string/236803#236803
template <typename Out>
void split(const std::string& s, char delim, Out result);
std::vector<std::string> split(const std::string& s, char delim);

// Taken from:
// https://stackoverflow.com/questions/1798112/removing-leading-and-trailing-spaces-from-a-string
std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v");
std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v");
std::string& trim(std::string& s, const char* t = " \t\n\r\f\v");

double toDouble(std::string& str);

/* 
 * Find point of intersection between a 1-D table polyline and a line segment 
 * using a linear search 
 * Returns:
 *  bool   = true  if the point of intersection is inside the 
 *                 table polyline and the line segment
 *         = false otherwise
 *  size_t = index of first point on intersecting table polyline segment 
 *  double = value of interpolation parameter t_p at the point of intersection
 *  double = value of interpolation parameter t_q at the point of intersection
 *  double  = x-coordinate of the point of intersection
 *  double  = y-coordinate of the point of intersection
 */
std::tuple<bool, std::size_t, double, double, double, double> 
  findIntersectionTableLinearSearch(const std::vector<double>& x_coords,
                                    const std::vector<double>& y_coords,
                                    double x_start_seg,
                                    double y_start_seg,
                                    double x_end_seg,
                                    double y_end_seg);

} // end namespace Vaango::Util
} // end namespace Vaango

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_UTIL_H
