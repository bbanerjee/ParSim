/*
 * The MIT License
 *
 * Copyright (c) 2017-2018 Parresia Research Limited, New Zealand
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

#ifndef GRADATION_H
#define GRADATION_H

#include <Core/Types/RealTypes.h>
#include <Core/Util/Utility.h>
#include <InputOutput/InputParameter.h>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <InputOutput/IOUtils.h>
#include <InputOutput/zenxml/xml.h>
#include <cstddef>
#include <vector>
#include <iostream>
#include <fstream>

namespace dem {

class Gradation
{
public:
  Gradation() 
    : sieveNum(0)
    , ratioBA(1)
    , ratioCA(1)
  {
  }

  Gradation(std::size_t sn, std::vector<REAL> v1, std::vector<REAL> v2, REAL ba,
            REAL ca)
    : sieveNum(sn)
    , percent(v1)
    , size(v2)
    , ratioBA(ba)
    , ratioCA(ca)
  {
  }

  inline void set(std::size_t sn, std::vector<REAL> v1, std::vector<REAL> v2,
                  REAL ba, REAL ca)
  {
    sieveNum = sn;
    percent = v1;
    size = v2;
    ratioBA = ba;
    ratioCA = ca;
  }

  /**
   * Initialize the particle size information from
   * CSV or XML files
   */
  void initializeFromCSVFile(std::ifstream& ifs);
  void initializeFromXMLFile(zen::XmlIn& ps);

  /**
   * Initialize the particle size information from
   * the global input parameters
   */
  void initializeFromInputParameters();

  /**
   * Output the gradation data in CSV or XMLformat
   */
  void outputCSV(std::ofstream& ofs) const;
  void outputXML(zen::XmlOut& xml) const;

  std::size_t getSieveNum() const { return sieveNum; }
  std::vector<REAL>& getPercent() { return percent; }
  const std::vector<REAL>& getPercent() const { return percent; }
  std::vector<REAL>& getSize() { return size; }
  const std::vector<REAL>& getSize() const { return size; }
  REAL getPtclRatioBA() const { return ratioBA; }
  REAL getPtclRatioCA() const { return ratioCA; }
  void setPtclRatioBA(REAL ba) { ratioBA = ba; }
  void setPtclRatioCA(REAL ca) { ratioCA = ca; }

  REAL getPtclMaxRadius() const { return (sieveNum > 0) ? size[0] : -1; }
  REAL getPtclMinRadius() const { return (sieveNum > 0) ? size[sieveNum - 1] * ratioCA : -1; }

  friend std::ostream& operator<<(std::ostream& out, const Gradation& grad)
  {
    zen::XmlDoc doc("Gradation");
    zen::XmlOut xml(doc);
    grad.outputXML(xml);
    out << zen::serialize(doc);
    return out;
  }

private:
  std::size_t sieveNum; // sieveNum == percent.size() == size.size()
  std::vector<REAL> percent;
  std::vector<REAL> size;
  REAL ratioBA; // ratio of radius b to radius a
  REAL ratioCA; // ratio of radius c to radius a

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& sieveNum;
    ar& percent;
    ar& size;
    ar& ratioBA;
    ar& ratioCA;
  }
};

} // namespace dem

#endif
