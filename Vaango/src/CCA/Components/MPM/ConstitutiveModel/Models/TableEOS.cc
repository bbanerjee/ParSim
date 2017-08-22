#include <CCA/Components/MPM/ConstitutiveModel/Models/TableEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/TableUtils.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>

using namespace Vaango;
using ProblemSpecP = Uintah::ProblemSpecP;
using ProblemSetupException = Uintah::ProblemSetupException;
using json = nlohmann::json;

TableEOS::TableEOS(ProblemSpecP& ps)
{
  std::string indep_vars;
  std::string dep_vars;
  ps->require("independent_variables", indep_vars);
  ps->require("dependent_variables", dep_vars);
  ps->require("filename", d_filename);

  std::vector<std::string> indepVarNames = parseVariableNames(indep_vars);
  std::vector<std::string> depVarNames = parseVariableNames(dep_vars);

  std::cout << "Independent:";
  int index = 0;
  for (const auto& name : indepVarNames) {
    d_indepVars.push_back(std::make_unique<IndependentVar>(name));
    std::cout << index << ":" << name << " ";
    index++;
  }
  std::cout << std::endl;

  std::cout << "Dependent:";
  index = 0;
  for (const auto& name : depVarNames) {
    d_depVars.push_back(std::make_unique<DependentVar>(name));
    std::cout << index << ":" << name << " ";
    index++;
  }
  std::cout << std::endl;
}

void 
TableEOS::addIndependentVariable(const std::string& varName)
{
  d_indepVars.push_back(std::make_unique<IndependentVar>(varName));
}

std::size_t  
TableEOS::addDependentVariable(const std::string& varName)
{
  // Check if dependent variable already exists
  // and just return index if true 
  for (auto ii = 0u; ii < d_depVars.size(); ii++) {
    if (d_depVars[ii]->name == varName) {
      return ii;
    }
  }

  // Name does not exist, add the new variable and return its index
  d_depVars.push_back(std::make_unique<DependentVar>(varName));
  return d_depVars.size()-1;  
}

void 
TableEOS::setup()
{
  readJSONTableFile(d_filename);
}

double 
TableEOS::interpolate(int index,
                      std::vector<double>& independents)
{
  return 0; 
}

json
TableEOS::loadJSON(std::stringstream& inputStream,
                   const std::string& tableFile)
{
  json doc;
  try {
    doc << inputStream;
  } catch (std::invalid_argument err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Cannot parse tabular input data file "
        << tableFile << "\n"
        << " Please check that the file is valid JSON using a linter";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return doc;
}

json
TableEOS::getContentsJSON(const json& doc,
                          const std::string& tableFile)
{
  json contents;
  try {
    contents = doc["Vaango_tabular_data"];
  } catch (std::out_of_range err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Cannot find the key \"Vaango_tabular_data\" in "
        << tableFile ;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return contents;
}

std::string
TableEOS::getTitleJSON(const json& contents,
                       const std::string& tableFile)
{
  std::string title;
  try {
    title = contents["Meta"].at("title").get<std::string>();
  } catch (std::out_of_range err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Could not find tabular EOS title in "
        << tableFile ;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return title;
}

json
TableEOS::getAllDataJSON(const json& contents,
                         const std::string& tableFile)
{
  json data;
  try {
    data = contents["Data"];
  } catch (std::exception err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Key \"Data\" not found in tabular EOS input file "
        << tableFile ;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return data;
}


std::vector<double>
TableEOS::getVectorJSON(const json& object,
                        const std::string key,
                        const std::string& tableFile)
{
  std::string str;
  try {
    str = object[key].get<std::string>();
  } catch (std::exception err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Dependent variable \""
        << key << "\" not found in tabular EOS input file "
        << tableFile ;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  std::vector<std::string> vecStr = Vaango::Util::split(str, ',');
  std::vector<double> vec;
  for (auto& val : vecStr) {
    vec.push_back(Vaango::Util::toDouble(val));
  }
  return vec;
}

void
TableEOS::readJSONTableFile(const std::string& tableFile)
{
  std::ifstream inputFile(tableFile);
  if (!inputFile) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Cannot read tabular input data file "
        << tableFile;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  std::stringstream inputStream;
  inputStream << inputFile.rdbuf();
  inputFile.close();

  json doc = loadJSON(inputStream, tableFile);
  json contents = getContentsJSON(doc, tableFile);
  std::string title = getTitleJSON(contents, tableFile);
  json data = getAllDataJSON(contents, tableFile);

  try {
    for (const auto& object : data) {
      std::string vecStr;
      try {
        vecStr = object[d_depVars[0]->name].get<std::string>();
      } catch (std::exception err) {
        std::ostringstream out;
        out << "**ERROR**"
            << " Dependent variable \""
            << d_depVars[0]->name << "\" not found in tabular EOS input file "
            << tableFile ;
        throw ProblemSetupException(out.str(), __FILE__, __LINE__);
      }
      json data0;
      try {
        data0 = data["Data"];
      } catch (std::exception err) {
        std::ostringstream out;
        out << "**ERROR**"
            << " Key \"Data\" not found in tabular EOS input file "
            << " for dependent variable \"" << d_depVars[0]->name << "\" in "
            << tableFile ;
        throw ProblemSetupException(out.str(), __FILE__, __LINE__);
      }
      for (const auto& object0 : data0) {
        std::string vecStr0;
        try {
          vecStr0 = object0[d_depVars[1]->name].get<std::string>();
        } catch (std::exception err) {
          std::ostringstream out;
          out << "**ERROR**"
              << " Dependent variable \""
              << d_depVars[1]->name << "\" not found in tabular EOS input file "
              << tableFile ;
          throw ProblemSetupException(out.str(), __FILE__, __LINE__);
        }
        json data1;
        try {
          data1 = data0["Data"];
        } catch (std::exception err) {
          std::ostringstream out;
          out << "**ERROR**"
              << " Key \"Data\" not found in tabular EOS input file "
              << " for dependent variable \"" << d_depVars[1]->name << "\" in "
              << tableFile ;
          throw ProblemSetupException(out.str(), __FILE__, __LINE__);
        }
        for (const auto& object1 : data1) {

          int index = 0;
          for (const auto& depVar : d_depVars) {
            if (index < 2) continue;

            std::string vecStr2;
            try {
              vecStr2 = object1[depVar->name].get<std::string>();
            } catch (std::exception err) {
              std::ostringstream out;
              out << "**ERROR**"
                  << " Dependent variable \""
                  << depVar->name << "\" not found in tabular EOS input file "
                  << tableFile ;
              throw ProblemSetupException(out.str(), __FILE__, __LINE__);
            }
            index++;
          }

          for (const auto& indepVar : d_indepVars) {
            std::string vecStr2;
            try {
              vecStr2 = object1[indepVar->name].get<std::string>();
            } catch (std::exception err) {
              std::ostringstream out;
              out << "**ERROR**"
                  << " Independent variable \""
                  << indepVar->name << "\" not found in tabular EOS input file "
                  << tableFile ;
              throw ProblemSetupException(out.str(), __FILE__, __LINE__);
            }
          }
        }
      }
    }
  } catch (std::exception err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Key \"Data\" not found in tabular EOS input file "
        << tableFile ;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
}

std::vector<std::string> 
TableEOS::parseVariableNames(const std::string& vars)
{
  std::vector<std::string> varNames = Vaango::Util::split(vars, ',');
  for (auto& name : varNames) {
    name = Vaango::Util::trim(name);
  }
  return varNames;
}

