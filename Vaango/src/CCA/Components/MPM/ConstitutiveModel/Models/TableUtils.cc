/*
 * The MIT License
 *
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/Models/TableUtils.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCondUtils.h>
#include <algorithm>
#include <sstream>

namespace Vaango {

namespace Util {

// Taken from:
// https://stackoverflow.com/questions/236129/most-elegant-way-to-split-a-string/236803#236803
template <typename Out>
void
split(const std::string& s, char delim, Out result)
{
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    *(result++) = item;
  }
}

std::vector<std::string>
split(const std::string& s, char delim)
{
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}

// trim from left
std::string&
ltrim(std::string& s, const char* t)
{
  s.erase(0, s.find_first_not_of(t));
  return s;
}

// trim from right
std::string&
rtrim(std::string& s, const char* t)
{
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

// trim from left & right
std::string&
trim(std::string& s, const char* t)
{
  return ltrim(rtrim(s, t), t);
}

// Convert into double
double
toDouble(std::string& str)
{
  // Remove illegal characters
  const std::string illegalChars = " [\\/:?\"<>|]";
  str.erase(remove_if(str.begin(), str.end(),
                      [&illegalChars](char cc) {
                        return (illegalChars.find(cc) != std::string::npos)
                                 ? true
                                 : false;
                      }),
            str.end());

  // Convert the string to double
  return std::stod(str);
}


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
                                  double y_end_seg)
{

  Uintah::Point P0(x_start_seg, y_start_seg, 0);
  Uintah::Point P1(x_end_seg, y_end_seg, 0);

  bool status = false;
  double t_p, t_q;
  Uintah::Point intersect;

  for (auto index = 0u; index < x_coords.size()-1; ++index) {
    Uintah::Point Q0(x_coords[index], y_coords[index], 0);
    Uintah::Point Q1(x_coords[index+1], y_coords[index+1], 0);

    std::tie(status, t_p, t_q, intersect) = 
      Vaango::Util::intersectionPoint(P0, P1, Q0, Q1);
    //std::cout << "P0 = " << P0 << " P1 = " << P1 
    //          << "Q0 = " << Q0 << " Q1 = " << Q1 << "\n";
    //std::cout << "status = " << status << " t_p = " << t_p << " t_q = " << t_q
    //          << " intersect = " << intersect << "\n";

    if (status) {
      return std::make_tuple(true, index, t_p, t_q, intersect.x(), intersect.y());
    } else {
      if (index == 0 && t_q < 0) {
        // Intersection is less than lowest value in table
        return std::make_tuple(true, index, t_p, t_q, intersect.x(), intersect.y());
      } else {
        if (index == x_coords.size()-2 && t_q > 1) {
          // Intersection is greater than largest value in table
          return std::make_tuple(true, index, t_p, t_q, intersect.x(), intersect.y());
        }
      }
    }
  }

  return std::make_tuple(false, x_coords.size(), t_p, t_q, intersect.x(), intersect.y());
}

} // End namespace Util

} // end namespace Vaango
