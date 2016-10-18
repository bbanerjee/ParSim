/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

/*
 * test_variant.cc
 *
 *  Created on: 18/10/2013
 *      Author: banerjee
 */

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <boost/variant.hpp>
#include <MPMDataTypes.h>

using namespace BrMPM;

typedef std::vector<double> VectorOfDouble;
typedef boost::variant<int, std::string, double, VectorOfDouble, std::vector<double> > variant_type;
typedef std::map<int, variant_type> map_type;
typedef std::map<std::string, variant_type> alt_map_type;

std::ostream&
operator<<(std::ostream& os, variant_type& var) {
  if (int* val = boost::get<int>(&var)) {
    os << " number " << *val;
  } else if (double* valDbl = boost::get<double>(&var)) {
    os << " number " << *valDbl;
  } else if (VectorOfDouble* valVec = boost::get<VectorOfDouble>(&var)) {
    os << " vector [";
    for (auto iter = valVec->begin(); iter != valVec->end(); ++iter) {
      os << *iter << " ";
    }
    os << "]";
  } else if (std::vector<double>* valvec = boost::get<std::vector<double> >(&var)) {
    os << " std::vector [";
    for (auto iter = valvec->begin(); iter != valvec->end(); ++iter) {
      os << *iter << " ";
    }
    os << "]";
  } else {
    std::string* valStr = boost::get<std::string>(&var);
    os << " string " << *valStr;
  }
  return os;
}

class PrintVisitor : public boost::static_visitor<std::string>
{
public:
  PrintVisitor(std::ostringstream& os) : d_os(os){}

  std::string operator()(int& val) {
    d_os << " visitor number " << val;
    return d_os.str();
  }
  std::string operator()(double& val) {
    d_os << " visitor number " << val;
    return d_os.str();
  }
  std::string operator()(VectorOfDouble& val) {
    d_os << " visitor vector [";
    for (auto iter = val.begin(); iter != val.end(); ++iter) {
      d_os << *iter << " ";
    }
    d_os << "]";
    return d_os.str();
  }
  std::string operator()(std::string& val) {
    d_os << " visitor string " << val;
    return d_os.str();
  }
private:
  std::ostringstream& d_os;
};

class LocalZeroVisitor : public boost::static_visitor<void>
{
public:
  void operator()(int& val) const {
    val = 0;
  }
  void operator()(double& val) const {
    val = 0.0;
  }
  void operator()(VectorOfDouble& val) const {
    for (auto iter = val.begin(); iter != val.end(); ++iter) {
      *iter = 0.0;
    }
  }
  void operator()(std::string& val) const {
    val = "0000";
  }
};

class GetVisitor : public boost::static_visitor<void>
{
public:
  template <typename T1>
  void operator()(T1& val) const
  {
  }
};

void add(alt_map_type& map, const std::string& label, variant_type& var)
{
  map[label] = var;
}

template<typename T>
void get(alt_map_type& map, const std::string& label, T& var)
{
  variant_type variant = map[label];
  var = boost::get<T>(variant);
}

template void get(alt_map_type& map, const std::string& label, std::vector<double>& var);
//{
//  variant_type variant = map[label];
//  var = boost::get<std::vector<double> >(variant);
//}


int main()
{
  // Set up initial data
  int data1 = 0;
  std::string data2 = "one";
  double data3 = 67.3;
  std::vector<double> data4;
  data4.emplace_back(1.5);
  data4.emplace_back(2.5);
  data4.emplace_back(3.5);

  // Create variants and fill map
  variant_type v1(data1);
  variant_type v2(data2);
  variant_type v3(data3);
  variant_type v4(data4);
  alt_map_type map1;
  map1["first"] = v1;
  map1["second"] = v2;
  map1["third"] = v3;
  map1["fourth"] = v4;
  double data5 = 7.91;
  variant_type v5(data5);
  add(map1, "fifth", v5);

  // print map
  for (auto iter = map1.begin(); iter != map1.end(); ++iter) {
    std::ostringstream os_map1;
    PrintVisitor pv_map1(os_map1);
    std::cout << "map1[" << iter->first << "]" << boost::apply_visitor(pv_map1, iter->second) << "\n";
  }
  std::cout << std::endl;

  // zero out the elements and print
  for (auto iter = map1.begin(); iter != map1.end(); ++iter) {
    boost::apply_visitor(LocalZeroVisitor(), iter->second);
    std::ostringstream os_map1;
    PrintVisitor pv_map1(os_map1);
    std::cout << "map1[" << iter->first << "]" << boost::apply_visitor(pv_map1, iter->second) << "\n";
  }
  std::cout << std::endl;

  VectorOfDouble var;
  //std::vector<double> var;
  get(map1, "fourth", var);
  std::cout << " print vector \"fourth\" [";
  for (auto iter = var.begin(); iter != var.end(); ++iter) {
    std::cout << *iter << " ";
  }
  std::cout << "]" << std::endl;


  map_type m;
  m[0] = variant_type(0);
  m[1] = variant_type("one");
  m[2] = variant_type(3.1);
  //m[0] = 0;
  //m[1] = "one";
  //m[2] = 3.1;

  std::vector<double> init;
  init.emplace_back(1.1);
  init.emplace_back(2.2);
  init.emplace_back(3.3);
  m[3] = init;

  VectorOfDouble vec(init);
  m[4] = vec;

  std::cout
  << "m[0] " << m[0] << "\n"
  << "m[1] " << m[1] << "\n"
  << "m[2] " << m[2] << "\n"
  << "m[3] " << m[3] << "\n"
  << "m[4] " << m[4] << "\n"
  << std::endl;

  std::ostringstream os;
  PrintVisitor pv(os);
  std::cout << "m[0] " << boost::apply_visitor(pv, m[0]) << "\n";
  std::ostringstream os1;
  PrintVisitor pv1(os1);
  std::cout << "m[1] " << boost::apply_visitor(pv1, m[1]) << "\n";
  std::ostringstream os2;
  PrintVisitor pv2(os2);
  std::cout << "m[2] " << boost::apply_visitor(pv2, m[2]) << "\n";
  std::ostringstream os3;
  PrintVisitor pv3(os3);
  std::cout << "m[3] " << boost::apply_visitor(pv3, m[3]) << "\n";
  std::ostringstream os4;
  PrintVisitor pv4(os4);
  std::cout << "m[4] " << boost::apply_visitor(pv4, m[4]) << "\n";
  std::cout << std::endl;

  boost::get<int>(m[0]) = 1;

  boost::get<std::string>(m[1]) = "new string";

  map_type::iterator iter = m.find(2);
  boost::get<double>(iter->second) = 5.2;

  init.push_back(5.5);
  boost::get<std::vector<double> >(m[3]) = init;

  vec.push_back(4.4);
  boost::get<VectorOfDouble>(m[4]) = vec;

  std::cout
  << "m[0] " << m[0] << "\n"
  << "m[1] " << m[1] << "\n"
  << "m[2] " << m[2] << "\n"
  << "m[3] " << m[3] << "\n"
  << "m[4] " << m[4] << "\n"
  << std::endl;

  std::cout << "m[0] " << boost::apply_visitor(pv, m[0]) << "\n";
  std::cout << "m[1] " << boost::apply_visitor(pv1, m[1]) << "\n";
  std::cout << "m[2] " << boost::apply_visitor(pv2, m[2]) << "\n";
  std::cout << "m[3] " << boost::apply_visitor(pv3, m[3]) << "\n";
  std::cout << "m[4] " << boost::apply_visitor(pv4, m[4]) << "\n";
  return 0;
}

