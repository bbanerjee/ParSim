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

#include <DiscreteElements/Gradation.h>
#include <Core/Const/Constants.h>
#include <iomanip>

using namespace dem;

void 
Gradation::initializeFromCSVFile(std::ifstream& ifs)
{
  std::size_t sieveNum;
  ifs >> sieveNum;
  std::vector<REAL> percent(sieveNum), size(sieveNum);
  for (auto i = 0u; i < sieveNum; ++i)
    ifs >> percent[i] >> size[i];
  REAL ratio_ba, ratio_ca;
  ifs >> ratio_ba >> ratio_ca;
  set(sieveNum, percent, size, ratio_ba, ratio_ca);
}

// **TODO** Add validity checks
void 
Gradation::initializeFromXMLFile(zen::XmlIn& ps)
{
  auto sieve_ps = ps["Sieves"];
  if (sieve_ps) {
    std::size_t numSieves;
    sieve_ps.attribute("number", numSieves);

    std::string percentPassingStr;
    sieve_ps["percent_passing"](percentPassingStr);
    std::vector<REAL> percentPassing = Ellip3D::IOUtil::convertStrArray<REAL>(percentPassingStr);

    std::string sizeStr;
    sieve_ps["size"](sizeStr);
    std::vector<REAL> size = Ellip3D::IOUtil::convertStrArray<REAL>(sizeStr);

    REAL ratio_ba, ratio_ca;
    sieve_ps["sieve_ratio"]["ratio_ba"](ratio_ba);
    sieve_ps["sieve_ratio"]["ratio_ca"](ratio_ca);

    set(numSieves, percentPassing, size, ratio_ba, ratio_ca);
  }
}

/**
 * Initialize the particle size information from
 * the global input parameters
 */
void 
Gradation::initializeFromInputParameters()
{
  auto sieveNum = util::getParam<std::size_t>("sieveNum");
  std::vector<REAL> percent(sieveNum), size(sieveNum);
  std::vector<std::pair<REAL, REAL>>& grada =
    InputParameter::get().gradation;
  assert(grada.size() == sieveNum);
  for (std::size_t i = 0; i < sieveNum; ++i) {
    percent[i] = grada[i].first;
    size[i] = grada[i].second;
  }
  REAL ratio_ba = util::getParam<REAL>("ratioBA");
  REAL ratio_ca = util::getParam<REAL>("ratioCA");

  set(sieveNum, percent, size, ratio_ba, ratio_ca);
}

/**
 * Output the gradation data in CSV format
 */
void 
Gradation::outputCSV(std::ofstream& ofs) const
{
  std::size_t sieveNum = getSieveNum();
  std::vector<REAL> percent = getPercent();
  std::vector<REAL> size = getSize();
  ofs << std::endl << std::setw(OWID) << sieveNum << std::endl;
  for (std::size_t i = 0; i < sieveNum; ++i)
    ofs << std::setw(OWID) << percent[i] << std::setw(OWID) << size[i]
        << std::endl;
  ofs << std::endl
      << std::setw(OWID) << getPtclRatioBA() << std::setw(OWID)
      << getPtclRatioCA() << std::endl;
}

/**
 * Output the gradation data in XML format
 */
void 
Gradation::outputXML(zen::XmlOut& xml) const
{
  std::size_t sieveNum = getSieveNum();
  std::vector<REAL> percent = getPercent();
  std::vector<REAL> size = getSize();

  auto& element = xml.ref();
  auto& child = element.addChild("Sieves");
  zen::XmlOut sieve_ps(child);
  sieve_ps.attribute("number", sieveNum);
  sieve_ps.attribute("compression", "none");
  sieve_ps.attribute("encoding", "ascii");

  std::ostringstream stream;
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.precision(dem::OPREC);
  for (const auto& value : percent) {
    stream << std::setw(OWID) << value;
  }
  sieve_ps["percent_passing"](stream.str());

  stream.str("");
  for (const auto& value : size) {
    stream << std::setw(OWID) << value;
  }
  sieve_ps["size"](stream.str());

  sieve_ps["sieve_ratio"]["ratio_ba"](getPtclRatioBA());
  sieve_ps["sieve_ratio"]["ratio_ca"](getPtclRatioCA());
}