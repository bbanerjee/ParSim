/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/XMLUtils.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <vector>

#include <libxml/tree.h>

using namespace Uintah;

ProblemSpec::ProblemSpec(const std::string& buffer)
  : d_documentNode(true)
{
  xmlDocPtr doc =
    xmlReadMemory(buffer.c_str(), buffer.length(), nullptr, nullptr, 0);
  d_node = xmlDocGetRootElement(doc);
}

ProblemSpecP
ProblemSpec::findBlock() const
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::findBlock()");
  const xmlNode* child = d_node->children;
  if (child != 0) {
    if (child->type == XML_TEXT_NODE) {
      child = child->next;
    }
  }
  if (child == nullptr) {
    return 0;
  } else {
    return scinew ProblemSpec(child, false);
  }
}

// Searchs through the given XML file 'fp' for a matching 'block'.
// Returns false if not found.  'fp' is updated to point to the line
// after the 'block'.
bool
ProblemSpec::findBlock(const std::string& name, FILE*& fp)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::findBlock(string,FILE)");

  while (true) {
    std::string line = UintahXML::getLine(fp);
    if (line == name) {
      return true;
    } else if (line == "") {
      return false;
    }
  }
}

ProblemSpecP
ProblemSpec::findBlock(const std::string& name) const
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::findBlock(string)");
  if (d_node == 0) {
    return 0;
  }
  const xmlNode* child = d_node->children;
  while (child != 0) {
    //    string child_name(to_char_ptr(child->name));
    std::string child_name((const char*)(child->name));
    if (name == child_name) {
      return scinew ProblemSpec(child, false);
    }
    child = child->next;
  }
  return 0;
}

//  This function recursively searches through the xml tree for xmlNodes
//  with the tagname and returns a vector of them.
std::vector<xmlNode*>
ProblemSpec::recursiveFind_xmlNodes(std::vector<xmlNode*> nodesFound,
                                    xmlNode* a_node,
                                    const string& tagname,
                                    const bool isTopLevelNode) const
{
  for (xmlNode* cur_node = a_node; cur_node; cur_node = cur_node->next) {

    string nodeName((const char*)(cur_node->name));
    bool isElementNode = (cur_node->type == XML_ELEMENT_NODE);
    bool isTagName     = (nodeName == tagname);

    if (isElementNode && isTagName) {
      nodesFound.push_back(cur_node);
    }

    xmlNode* child = cur_node->children;
    if (child != 0) {
      if (child->type == XML_ELEMENT_NODE) {
        nodesFound = recursiveFind_xmlNodes(nodesFound, child, tagname, false);
      }
    }

    // Avoid looping over the siblings of top level node
    if (isTopLevelNode) {
      break;
    }
  }
  return nodesFound;
}

//______________________________________________________________________
//  This function recursively searches through the xml tree for xmlNodes
//  with the tagname and returns a vector ProblemSpecs.
std::vector<ProblemSpecP>
ProblemSpec::findBlocksRecursive(const string& tagname) const
{
  std::vector<ProblemSpecP> results;
  if (d_node == 0) {
    return results;
  }

  std::vector<xmlNode*> xmlNodeResults;
  xmlNode* cur_node = d_node;

  xmlNodeResults =
    recursiveFind_xmlNodes(xmlNodeResults, cur_node, tagname, true);

  for (size_t i = 0; i < xmlNodeResults.size(); i++) {
    ProblemSpecP ps = scinew ProblemSpec(xmlNodeResults[i], false);
    results.push_back(ps);
  }
  return results;
}

ProblemSpecP
ProblemSpec::findBlockWithAttribute(const std::string& name,
                                    const std::string& attribute) const
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::findBlockWithAttribute(string,string)");

  for (ProblemSpecP ps = this->findBlock(name); ps != 0;
       ps              = ps->findNextBlock(name)) {

    std::string attr = "";
    ps->getAttribute(attribute, attr);
    if (attr.length() > 0) {
      return ps;
    } else {
      continue;
    }
  }

  return 0;
}

//  Finds:  <Block attribute = "value">
ProblemSpecP
ProblemSpec::findBlockWithAttributeValue(const string& name,
                                         const string& attribute,
                                         const string& value) const
{
  for (ProblemSpecP ps = this->findBlock(name); ps != nullptr;
       ps              = ps->findNextBlock(name)) {

    string attr = "";
    ps->getAttribute(attribute, attr);

    if (attr == value) {
      return ps;
    } else {
      continue;
    }
  }

  return nullptr;
}

ProblemSpecP
ProblemSpec::findBlockWithOutAttribute(const std::string& name) const

{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::findBlockWithOutAttribute(string)");

  for (ProblemSpecP ps = this->findBlock(name); ps != 0;
       ps              = ps->findNextBlock(name)) {

    std::map<std::string, std::string> attributes;
    ps->getAttributes(attributes);
    if (attributes.empty()) {
      return ps;
    } else {
      continue;
    }
  }

  return 0;
}

ProblemSpecP
ProblemSpec::findNextBlock() const
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::findNextBlock()");
  const xmlNode* found_node = d_node->next;

  if (found_node != 0) {
    if (found_node->type == XML_TEXT_NODE) {
      found_node = found_node->next;
    }
  }

  if (found_node == nullptr) {
    return 0;
  } else {
    return scinew ProblemSpec(found_node, false);
  }
}

ProblemSpecP
ProblemSpec::findNextBlock(const std::string& name) const
{
  MALLOC_TRACE_TAG_SCOPE("findNextBlock(string)");
  // Iterate through all of the child nodes of the next node
  // until one is found that has this name

  const xmlNode* found_node = d_node->next;

  while (found_node != 0) {
    //    string c_name(to_char_ptr(found_node->name));
    std::string c_name((const char*)(found_node->name));
    if (c_name == name) {
      break;
    }

    found_node = found_node->next;
  }
  if (found_node == nullptr) {
    return 0;
  } else {
    return scinew ProblemSpec(found_node, false);
  }
}

ProblemSpecP
ProblemSpec::findTextBlock()
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::findTextBlock()");
  for (xmlNode* child = d_node->children; child != 0; child = child->next) {
    if (child->type == XML_TEXT_NODE) {
      return scinew ProblemSpec(child, false);
    }
  }
  return nullptr;
}

std::string
ProblemSpec::getNodeName() const
{
  //  return string(to_char_ptr(d_node->name));
  return std::string((const char*)(d_node->name));
}

short
ProblemSpec::getNodeType()
{
  return d_node->type;
}

ProblemSpecP
ProblemSpec::importNode(ProblemSpecP src, bool deep)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::ImportNode()");
  xmlNode* d = xmlDocCopyNode(src->d_node, d_node->doc, deep ? 1 : 0);
  if (d) {
    return scinew ProblemSpec(d, false);
  } else {
    return 0;
  }
}

void
ProblemSpec::addComment(std::string comment)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::addComment()");
  xmlNodePtr commentNode = xmlNewComment(BAD_CAST comment.c_str());
  xmlAddChild(d_node, commentNode);
}

ProblemSpecP
ProblemSpec::makeComment(std::string comment)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::makeComment()");
  xmlNodePtr commentNode = xmlNewComment(BAD_CAST comment.c_str());
  return scinew ProblemSpec(commentNode, false);
}

void
ProblemSpec::replaceChild(ProblemSpecP toreplace, ProblemSpecP replaced)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::replaceChild()");
  xmlNode* d = xmlReplaceNode(toreplace->d_node, replaced->d_node);

  if (d) {
    xmlFreeNode(d);
  }
}

void
ProblemSpec::removeChild(ProblemSpecP child)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::removeChild()");
  xmlUnlinkNode(child->getNode());
  xmlFreeNode(child->getNode());
}

ProblemSpecP
ProblemSpec::get(const std::string& name, double& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  ProblemSpecP ps;

  std::string stringValue;
  ps = get(name, stringValue);
  if (ps == 0) {
    return ps;
  } else {
    UintahXML::validateType(stringValue, UintahXML::FLOAT_TYPE);
    std::istringstream ss(stringValue);
    ss >> value;
    if (!ss) {
      ps = 0;
      //      std::cout << "WARNING: ProblemSpec.cc: get(%s, double):
      //       std::stringstream failed..." << name << std::endl;
    }
  }

  return ps;
}

ProblemSpecP
ProblemSpec::get(const std::string& name, unsigned int& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  ProblemSpecP ps;

  std::string stringValue;
  ps = get(name, stringValue);
  if (ps == 0) {
    return ps;
  } else {
    UintahXML::validateType(stringValue, UintahXML::INT_TYPE);
    std::istringstream ss(stringValue);
    ss >> value;
    if (!ss) {
      printf("WARNING: ProblemSpec.cc: get(%s, uint):  std::stringstream "
             "failed...\n",
             name.c_str());
      ps = 0;
    }
  }

  return ps;
}

ProblemSpecP
ProblemSpec::get(const std::string& name, int& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  ProblemSpecP ps;

  std::string stringValue;
  ps = get(name, stringValue);
  if (ps == 0) {
    return ps;
  } else {
    UintahXML::validateType(stringValue, UintahXML::INT_TYPE);
    std::istringstream ss(stringValue);
    ss >> value;
    if (!ss) {
      printf(
        "WARNING: ProblemSpec.cc: get(%s, int):  std::stringstream failed...\n",
        name.c_str());
      ps = 0;
    }
  }

  return ps;
}

ProblemSpecP
ProblemSpec::get(const std::string& name, long& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  ProblemSpecP ps;

  std::string stringValue;
  ps = get(name, stringValue);
  if (ps == 0) {
    return ps;
  } else {
    UintahXML::validateType(stringValue, UintahXML::INT_TYPE);
    std::istringstream ss(stringValue);
    ss >> value;
    if (!ss) {
      printf("WARNING: ProblemSpec.cc: get(%s, long):  std::stringstream "
             "failed...\n",
             name.c_str());
      ps = 0;
    }
  }

  return ps;
}

ProblemSpecP
ProblemSpec::get(const std::string& name, bool& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  ProblemSpecP ps;

  std::string stringValue;
  ps = get(name, stringValue);
  if (ps == 0) {
    return ps;
  } else {
    // Slurp up any spaces that were put in before or after the cmp string.
    std::istringstream result_stream(stringValue);
    std::string nospace_cmp;
    result_stream >> nospace_cmp;

    if (!result_stream) {
      printf("WARNING: ProblemSpec.cc: get(%s, bool):  std::stringstream "
             "failed...\n",
             name.c_str());
    }

    if (nospace_cmp == "false") {
      value = false;
    } else if (nospace_cmp == "true") {
      value = true;
    } else {
      std::string error = name + " Must be either true or false";
      throw ProblemSetupException(error, __FILE__, __LINE__);
    }
  }
  return ps;
}

ProblemSpecP
ProblemSpec::get(const std::string& name, std::string& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  // the other gets will call this one to get the string...
  ProblemSpecP ps   = this;
  ProblemSpecP node = findBlock(name);
  if (node == 0) {
    ps = 0;
    return ps;
  } else { // eliminate spaces
    value = node->getNodeValue();

    // elminate spaces from string

    std::stringstream in_stream(value);
    std::vector<std::string> vs;
    copy(std::istream_iterator<std::string>(in_stream),
         std::istream_iterator<std::string>(),
         back_inserter(vs));
    std::string out_string;
    for (std::vector<std::string>::const_iterator it = vs.begin();
         it != vs.end();
         ++it) {
      out_string += *it + ' ';
    }

    if (out_string.length() > 0) {
      // if user accidentally leaves out value, this will crash with an ugly
      // exception
      std::string::iterator begin = out_string.end() - 1;
      std::string::iterator end   = out_string.end();
      out_string.erase(begin, end);
    }
    value = out_string;
  }
  return ps;
}

ProblemSpecP
ProblemSpec::get(const std::string& name, Point& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  Vector v;
  ProblemSpecP ps = get(name, v);
  value           = Point(v);
  return ps;
}

ProblemSpecP
ProblemSpec::get(const std::string& name, Vector& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  ProblemSpecP ps;

  std::string stringValue;
  ps = get(name, stringValue);
  if (ps == 0) {
    return ps;
  } else {
    // Parse out the [num,num,num]
    // Now pull apart the stringValue
    value = Vector::fromString(stringValue);
  }

  return ps;
}

ProblemSpec::InputType
ProblemSpec::getInputType(const std::string& stringValue)
{
  std::string validChars(" +-.0123456789eE");
  string::size_type pos = stringValue.find_first_not_of(validChars);
  if (pos != std::string::npos) {
    // we either have a string or a vector
    if (stringValue.find_first_of("[") == 0) {
      // this is most likely a vector vector
      return ProblemSpec::VECTOR_TYPE;
    } else {
      // we have a string
      return ProblemSpec::STRING_TYPE;
    }
  } else {
    // otherwise we have a number
    return ProblemSpec::NUMBER_TYPE;
  }
  return ProblemSpec::UNKNOWN_TYPE;
}

// value should probably be empty before calling this...
ProblemSpecP
ProblemSpec::get(const std::string& name, std::vector<double>& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  std::vector<std::string> string_values;
  if (!this->get(name, string_values)) {
    return 0;
  }

  for (std::vector<std::string>::const_iterator vit(string_values.begin());
       vit != string_values.end();
       vit++) {
    const std::string v(*vit);

    UintahXML::validateType(v, UintahXML::FLOAT_TYPE);
    value.push_back(atof(v.c_str()));
  }

  return this;
}

// value should probably be empty before calling this...
ProblemSpecP
ProblemSpec::get(const std::string& name,
                 std::vector<double>& value,
                 const int nItems)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  std::vector<std::string> string_values;
  if (!this->get(name, string_values, nItems)) {
    return 0;
  }

  for (std::vector<std::string>::const_iterator vit(string_values.begin());
       vit != string_values.end();
       vit++) {
    const std::string v(*vit);

    UintahXML::validateType(v, UintahXML::FLOAT_TYPE);
    value.push_back(atof(v.c_str()));
  }

  return this;
}

// value should probably be empty before calling this...
ProblemSpecP
ProblemSpec::get(const std::string& name, std::vector<int>& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  std::vector<std::string> string_values;
  if (!this->get(name, string_values)) {
    return 0;
  }

  for (std::vector<std::string>::const_iterator vit(string_values.begin());
       vit != string_values.end();
       vit++) {
    const std::string v(*vit);

    UintahXML::validateType(v, UintahXML::INT_TYPE);
    value.push_back(atoi(v.c_str()));
  }

  return this;
}

// value should probably be empty before calling this...
ProblemSpecP
ProblemSpec::get(const std::string& name, std::vector<std::string>& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  ProblemSpecP ps;

  std::string stringValue;
  ps = get(name, stringValue);
  if (ps == 0) {
    return ps;
  } else {
    std::istringstream in(stringValue);
    char c, next;
    std::string result;
    while (!in.eof()) {
      in >> c;
      if (c == '[' || c == ',' || c == ' ' || c == ']') {
        continue;
      }
      next = in.peek();
      result += c;
      if (next == ',' || next == ' ' || next == ']' || in.eof()) {
        // push next string onto stack
        value.push_back(result);
        result.erase();
      }
    }
  }
  return ps;
}

// value should probably be empty before calling this...
ProblemSpecP
ProblemSpec::get(const std::string& name,
                 std::vector<std::string>& value,
                 const int nItems)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  ProblemSpecP ps;

  std::string stringValue;
  ps = get(name, stringValue);
  if (ps == 0) {
    return ps;
  } else {
    std::istringstream in(stringValue);
    char c, next;
    std::string result;
    int counter = 0;
    while (!in.eof() && counter < nItems) {
      in >> c;
      if (c == '[' || c == ',' || c == ' ' || c == ']') {
        continue;
      }
      next = in.peek();
      result += c;
      if (next == ',' || next == ' ' || next == ']' || in.eof()) {
        // push next string onto stack
        value.push_back(result);
        result.erase();
        counter++;
      }
    }
  }
  return ps;
}

ProblemSpecP
ProblemSpec::get(const std::string& name, IntVector& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  ProblemSpecP ps;

  std::string stringValue;
  ps = get(name, stringValue);
  if (ps != 0) {
    value = IntVector::fromString(stringValue);
  }

  return ps;
}

ProblemSpecP
ProblemSpec::get(const std::string& name, std::vector<IntVector>& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  ProblemSpecP ps;

  std::string stringValue;
  ps = get(name, stringValue);
  if (ps == 0) {
    return ps;
  } else {
    std::istringstream in(stringValue);
    char c;
    bool first_bracket = false;
    bool inner_bracket = false;
    std::string result;
    // what we're going to do is look for the first [ then pass it.
    // then if we find another [, make a string out of it until we see ],
    // then pass that into parseIntVector, and repeat.
    while (!in.eof()) {
      in >> c;
      if (c == ' ' || (c == ',' && !inner_bracket)) {
        continue;
      }
      if (c == '[') {
        if (!first_bracket) {
          first_bracket = true;
          continue;
        } else {
          inner_bracket = true;
        }
      } else if (c == ']') {
        if (inner_bracket) {
          // parse the string for an IntVector
          IntVector val;
          result += c;
          // it should be [num,num,num] by now
          val = IntVector::fromString(result);
          value.push_back(val);
          result.erase();
          inner_bracket = false;
          continue;
        } else {
          break; // end parsing on outer ]
        }
      }
      // add the char to the string
      result += c;
    } // end while (!in.eof())
  }

  return ps;
}

bool
ProblemSpec::get(int& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  std::string stringValue;
  if (!get(stringValue)) {
    return false;
  }
  value = atoi(stringValue.c_str());
  return true;
}

bool
ProblemSpec::get(long& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  std::string stringValue;
  if (!get(stringValue)) {
    return false;
  }
  value = atoi(stringValue.c_str());
  return true;
}

bool
ProblemSpec::get(double& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  std::string stringValue;
  if (!get(stringValue)) {
    return false;
  }
  value = atof(stringValue.c_str());
  return true;
}

bool
ProblemSpec::get(std::string& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  std::string tmp = getNodeValue();
  if (tmp == "") {
    return false;
  }
  std::istringstream tmp_str(tmp);
  std::string w;
  while (tmp_str >> w) {
    value += w;
  }
  return true;
}

bool
ProblemSpec::get(Vector& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::get()");
  std::string stringValue;
  if (!get(stringValue)) {
    return false;
  }
  // Now pull apart the stringValue
  std::string::size_type i1 = stringValue.find("[");
  std::string::size_type i2 = stringValue.find_first_of(",");
  std::string::size_type i3 = stringValue.find_last_of(",");
  std::string::size_type i4 = stringValue.find("]");

  std::string x_val(stringValue, i1 + 1, i2 - i1 - 1);
  std::string y_val(stringValue, i2 + 1, i3 - i2 - 1);
  std::string z_val(stringValue, i3 + 1, i4 - i3 - 1);

  UintahXML::validateType(x_val, UintahXML::FLOAT_TYPE);
  UintahXML::validateType(y_val, UintahXML::FLOAT_TYPE);
  UintahXML::validateType(z_val, UintahXML::FLOAT_TYPE);

  value.x(atof(x_val.c_str()));
  value.y(atof(y_val.c_str()));
  value.z(atof(z_val.c_str()));
  return true;
}

ProblemSpecP
ProblemSpec::getWithDefault(const std::string& name,
                            double& value,
                            double defaultVal)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::getWithDefault()");
  ProblemSpecP ps = get(name, value);
  if (ps == 0) {

    // create xmlNode to add to the tree
    appendElement(name.c_str(), defaultVal);

    // set default values
    ps    = this;
    value = defaultVal;
  }

  return ps;
}

ProblemSpecP
ProblemSpec::getWithDefault(const std::string& name, int& value, int defaultVal)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::getWithDefault()");
  ProblemSpecP ps = get(name, value);
  if (ps == 0) {

    // create xmlNode to add to the tree
    appendElement(name.c_str(), defaultVal);

    // set default values
    ps    = this;
    value = defaultVal;
  }

  return ps;
}
ProblemSpecP
ProblemSpec::getWithDefault(const std::string& name,
                            bool& value,
                            bool defaultVal)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::getWithDefault()");
  ProblemSpecP ps = get(name, value);
  if (ps == 0) {

    // create xmlNode to add to the tree
    appendElement(name.c_str(), defaultVal);

    // set default values
    ps    = this;
    value = defaultVal;
  }

  return ps;
}
ProblemSpecP
ProblemSpec::getWithDefault(const std::string& name,
                            std::string& value,
                            const std::string& defaultVal)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::getWithDefault()");
  ProblemSpecP ps = get(name, value);
  if (ps == 0) {

    // create xmlNode to add to the tree
    appendElement(name.c_str(), defaultVal);

    // set default values
    ps    = this;
    value = defaultVal;
  }
  return ps;
}
ProblemSpecP
ProblemSpec::getWithDefault(const std::string& name,
                            IntVector& value,
                            const IntVector& defaultVal)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::getWithDefault()");
  ProblemSpecP ps = get(name, value);
  if (ps == 0) {

    // create xmlNode to add to the tree
    appendElement(name.c_str(), defaultVal);

    // set default values
    ps    = this;
    value = defaultVal;
  }

  return ps;
}

ProblemSpecP
ProblemSpec::getWithDefault(const std::string& name,
                            Vector& value,
                            const Vector& defaultVal)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::getWithDefault()");
  ProblemSpecP ps = get(name, value);
  if (ps == 0) {

    // create xmlNode to add to the tree
    appendElement(name.c_str(), defaultVal);

    // set default values
    ps    = this;
    value = defaultVal;
  }

  return ps;
}

ProblemSpecP
ProblemSpec::getWithDefault(const std::string& name,
                            Point& value,
                            const Point& defaultVal)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::getWithDefault()");
  ProblemSpecP ps = get(name, value);
  if (ps == 0) {

    // create xmlNode to add to the tree
    appendElement(name.c_str(), defaultVal);

    // set default values
    ps    = this;
    value = defaultVal;
  }

  return ps;
}

ProblemSpecP
ProblemSpec::getWithDefault(const std::string& name,
                            std::vector<double>& value,
                            const std::vector<double>& defaultVal)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::getWithDefault()");
  value.clear();
  ProblemSpecP ps = get(name, value);
  if (ps == 0) {

    // create xmlNode to add to the tree
    appendElement(name.c_str(), defaultVal);

    // set default values
    ps = this;

    value.clear();
    int size = static_cast<int>(defaultVal.size());
    for (int i = 0; i < size; i++) {
      value.push_back(defaultVal[i]);
    }
  }

  return ps;
}

ProblemSpecP
ProblemSpec::getWithDefault(const std::string& name,
                            std::vector<int>& value,
                            const std::vector<int>& defaultVal)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::getWithDefault()");
  value.clear();
  ProblemSpecP ps = get(name, value);
  if (ps == 0) {

    // add xmlNode to the tree
    appendElement(name.c_str(), defaultVal);
    // set default values
    ps = this;
    value.clear();
    int size = static_cast<int>(defaultVal.size());
    for (int i = 0; i < size; i++) {
      value.push_back(defaultVal[i]);
    }
  }

  return ps;
}

ProblemSpecP
ProblemSpec::appendElement(const char* name, const std::string& value)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::appendElement()");
  xmlNode* newnode =
    xmlNewChild(d_node, 0, BAD_CAST name, BAD_CAST value.c_str());
  return scinew ProblemSpec(newnode, false);
}

// basically to make sure correct overloaded function is called
ProblemSpecP
ProblemSpec::appendElement(const char* name, const char* value)
{
  return appendElement(name, std::string(value));
}

ProblemSpecP
ProblemSpec::appendElement(const char* name, int value)
{
  std::ostringstream val;
  val << value;
  return appendElement(name, val.str());
}

ProblemSpecP
ProblemSpec::appendElement(const char* name, unsigned int value)
{
  std::ostringstream val;
  val << value;
  return appendElement(name, val.str());
}

ProblemSpecP
ProblemSpec::appendElement(const char* name, long value)
{
  std::ostringstream val;
  val << value;
  return appendElement(name, val.str());
}

ProblemSpecP
ProblemSpec::appendElement(const char* name, const IntVector& value)
{
  std::ostringstream val;
  val << '[' << value.x() << ", " << value.y() << ", " << value.z() << ']';
  return appendElement(name, val.str());
}

ProblemSpecP
ProblemSpec::appendElement(const char* name, const Point& value)
{

  std::ostringstream val;
  val << '[' << std::setprecision(17) << value.x() << ", "
      << std::setprecision(17) << value.y() << ", " << std::setprecision(17)
      << value.z() << ']';
  return appendElement(name, val.str());
}

ProblemSpecP
ProblemSpec::appendElement(const char* name, const Vector& value)
{
  std::ostringstream val;
  val << '[' << std::setprecision(17) << value.x() << ", "
      << std::setprecision(17) << value.y() << ", " << std::setprecision(17)
      << value.z() << ']';
  return appendElement(name, val.str());
}

ProblemSpecP
ProblemSpec::appendElement(const char* name, double value)
{
  std::ostringstream val;
  val << std::setprecision(17) << value;
  return appendElement(name, val.str());
}

ProblemSpecP
ProblemSpec::appendElement(const char* name, const std::vector<double>& value)
{
  std::ostringstream val;
  val << '[';
  for (unsigned int i = 0; i < value.size(); i++) {
    val << std::setprecision(17) << value[i];
    if (i != value.size() - 1) {
      val << ',';
    }
  }
  val << ']';
  return appendElement(name, val.str());
}

ProblemSpecP
ProblemSpec::appendElement(const char* name, const std::vector<int>& value)
{
  std::ostringstream val;
  val << '[';
  for (unsigned int i = 0; i < value.size(); i++) {
    val << std::setprecision(17) << value[i];
    if (i != value.size() - 1) {
      val << ',';
    }
  }
  val << ']';
  return appendElement(name, val.str());
}

ProblemSpecP
ProblemSpec::appendElement(const char* name,
                           const std::vector<std::string>& value)
{
  std::ostringstream val;
  val << '[';
  for (unsigned int i = 0; i < value.size(); i++) {
    val << value[i];
    if (i != value.size() - 1) {
      val << ',';
    }
  }
  val << ']';
  return appendElement(name, val.str());
}

ProblemSpecP
ProblemSpec::appendElement(const char* name, bool value)
{
  if (value) {
    return appendElement(name, std::string("true"));
  } else {
    return appendElement(name, std::string("false"));
  }
}

void
ProblemSpec::require(const std::string& name, double& value)
{
  // Check if the prob_spec is nullptr
  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

void
ProblemSpec::require(const std::string& name, int& value)
{
  // Check if the prob_spec is nullptr
  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

void
ProblemSpec::require(const std::string& name, unsigned int& value)
{
  // Check if the prob_spec is nullptr
  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

void
ProblemSpec::require(const std::string& name, long& value)
{
  // Check if the prob_spec is nullptr
  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

void
ProblemSpec::require(const std::string& name, bool& value)
{
  // Check if the prob_spec is nullptr
  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

void
ProblemSpec::require(const std::string& name, std::string& value)
{
  // Check if the prob_spec is nullptr
  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

void
ProblemSpec::require(const std::string& name, Vector& value)
{
  // Check if the prob_spec is nullptr
  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

void
ProblemSpec::require(const std::string& name, std::vector<double>& value)
{

  // Check if the prob_spec is nullptr

  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

void
ProblemSpec::require(const std::string& name, std::vector<std::string>& value)
{

  // Check if the prob_spec is nullptr

  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

void
ProblemSpec::require(const std::string& name, std::vector<int>& value)
{

  // Check if the prob_spec is nullptr

  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

void
ProblemSpec::require(const std::string& name, std::vector<IntVector>& value)
{
  // Check if the prob_spec is nullptr
  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

void
ProblemSpec::require(const std::string& name, IntVector& value)
{
  // Check if the prob_spec is nullptr
  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

void
ProblemSpec::require(const std::string& name, Point& value)
{
  // Check if the prob_spec is nullptr
  if (!this->get(name, value)) {
    throw ParameterNotFound(name, __FILE__, __LINE__);
  }
}

bool
ProblemSpec::findAttribute(const std::string& attribute) const
{
  std::map<std::string, std::string> attributes;
  getAttributes(attributes);

  std::map<std::string, std::string>::iterator iter =
    attributes.find(attribute);

  if (iter != attributes.end()) {
    return true;
  } else {
    return false;
  }
}

void
ProblemSpec::getAttributes(std::map<std::string, std::string>& attributes) const
{
  attributes.clear();

  xmlAttr* attr = d_node->properties;

  for (; attr != 0; attr = attr->next) {
    if (attr->type == XML_ATTRIBUTE_NODE) {
      attributes[(const char*)(attr->name)] =
        (const char*)(attr->children->content);
    }
  }
}

bool
ProblemSpec::getAttribute(const std::string& attribute,
                          std::string& result) const
{

  std::map<std::string, std::string> attributes;
  getAttributes(attributes);

  std::map<std::string, std::string>::iterator iter =
    attributes.find(attribute);

  if (iter != attributes.end()) {
    result = iter->second;
    return true;
  } else {
    return false;
  }
}

//______________________________________________________________________
//
bool
ProblemSpec::getAttribute(const std::string& name,
                          std::vector<double>& value) const
{
  std::vector<std::string> stringValues;
  if (!getAttribute(name, stringValues)) {
    return false;
  }
  for (std::vector<std::string>::const_iterator vit(stringValues.begin());
       vit != stringValues.end();
       vit++) {
    const std::string v(*vit);
    UintahXML::validateType(v, UintahXML::FLOAT_TYPE);
    value.push_back(atof(v.c_str()));
  }
  return true;
}

//______________________________________________________________________
//
bool
ProblemSpec::getAttribute(const std::string& name, Vector& value) const
{
  std::string stringValue;
  if (!getAttribute(name, stringValue)) {
    return false;
  }
  // Parse out the [num,num,num]
  // Now pull apart the stringValue
  std::string::size_type i1 = stringValue.find("[");
  std::string::size_type i2 = stringValue.find_first_of(",");
  std::string::size_type i3 = stringValue.find_last_of(",");
  std::string::size_type i4 = stringValue.find("]");

  std::string x_val(stringValue, i1 + 1, i2 - i1 - 1);
  std::string y_val(stringValue, i2 + 1, i3 - i2 - 1);
  std::string z_val(stringValue, i3 + 1, i4 - i3 - 1);

  UintahXML::validateType(x_val, UintahXML::FLOAT_TYPE);
  UintahXML::validateType(y_val, UintahXML::FLOAT_TYPE);
  UintahXML::validateType(z_val, UintahXML::FLOAT_TYPE);

  value.x(atof(x_val.c_str()));
  value.y(atof(y_val.c_str()));
  value.z(atof(z_val.c_str()));

  return true;
}

bool
ProblemSpec::getAttribute(const std::string& name, double& value) const
{
  std::string stringValue;
  if (!getAttribute(name, stringValue)) {
    return false;
  }
  UintahXML::validateType(stringValue, UintahXML::FLOAT_TYPE);
  std::istringstream ss(stringValue);
  ss >> value;
  if (!ss) {
    printf(
      "WARNING: ProblemSpec.cc: getAttribute(%s, double):  std::stringstream "
      "failed...\n",
      name.c_str());
  }

  return true;
}

//______________________________________________________________________
//
bool
ProblemSpec::getAttribute(const std::string& name, int& value) const
{
  std::string stringValue;
  if (!getAttribute(name, stringValue)) {
    return false;
  }
  UintahXML::validateType(stringValue, UintahXML::INT_TYPE);
  std::istringstream ss(stringValue);
  ss >> value;
  if (!ss) {
    printf("WARNING: ProblemSpec.cc: getAttribute(%s, int):  std::stringstream "
           "failed...\n",
           name.c_str());
  }

  return true;
}

//______________________________________________________________________
//
bool
ProblemSpec::getAttribute(const std::string& name, bool& value) const
{
  std::string stringValue;
  if (!getAttribute(name, stringValue)) {
    return false;
  }
  // remove any spaces that were put in before or after the cmp string.
  std::istringstream result_stream(stringValue);
  std::string nospace_cmp;
  result_stream >> nospace_cmp;
  if (!result_stream) {
    printf(
      "WARNING: ProblemSpec.cc: get(%s, bool):  std::stringstream failed...\n",
      name.c_str());
  }

  if (nospace_cmp == "false") {
    value = false;
  } else if (nospace_cmp == "true") {
    value = true;
  } else {
    std::string error = name + " Must be either true or false";
    throw ProblemSetupException(error, __FILE__, __LINE__);
  }
  return true;
}

bool
ProblemSpec::getAttribute(const std::string& attribute,
                          std::vector<std::string>& result) const
{
  std::map<std::string, std::string> attributes;
  getAttributes(attributes);

  std::map<std::string, std::string>::iterator iter =
    attributes.find(attribute);

  if (iter != attributes.end()) {
    std::string attributeName = iter->second;
    std::stringstream ss(attributeName);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    result.assign(begin, end);
    // std::std::vector<std::string> vstrings(begin, end);
    // result = iter->second;
    return true;
  } else {
    return false;
  }
}

void
ProblemSpec::setAttribute(const std::string& name, const std::string& value)
{
  xmlNewProp(d_node, BAD_CAST name.c_str(), BAD_CAST value.c_str());
}

//______________________________________________________________________
//
void
ProblemSpec::removeAttribute(const std::string& attrName)
{
  xmlAttr* attr = xmlHasProp(d_node, BAD_CAST attrName.c_str());
  if (attr) {
    xmlRemoveProp(attr);
  }
}
//______________________________________________________________________
//
void
ProblemSpec::replaceAttributeValue(const std::string& attrName,
                                   const std::string& newValue)
{
  removeAttribute(attrName);
  setAttribute(attrName, newValue);
}

ProblemSpecP
ProblemSpec::getFirstChild()
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::getFirstChild()");
  xmlNode* d = d_node->children;
  if (d) {
    return scinew ProblemSpec(d, false);
  } else {
    return 0;
  }
}

ProblemSpecP
ProblemSpec::getNextSibling()
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::getNextSibling()");
  xmlNode* d = d_node->next;
  if (d) {
    return scinew ProblemSpec(d, false);
  } else {
    return 0;
  }
}

ProblemSpecP
ProblemSpec::getParent()
{
  xmlNode* d = d_node->parent;

  if (d) {
    return scinew ProblemSpec(d, false);
  } else {
    return nullptr;
  }
}

std::string
ProblemSpec::getNodeValue()
{
  std::string ret;
  for (xmlNode* child = d_node->children; child != 0; child = child->next) {
    if (child->type == XML_TEXT_NODE) {
      ret = (const char*)(child->content);
      break;
    }
  }
  return ret;
}

// append element with associated string
ProblemSpecP
ProblemSpec::appendChild(const char* str)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::appendChild()");
  xmlNode* elt = xmlNewChild(d_node, 0, BAD_CAST str, 0);

  return scinew ProblemSpec(elt, false);
}

void
ProblemSpec::appendChild(ProblemSpecP pspec)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::appendChild()");
  xmlAddChild(d_node, pspec->d_node);
}

// filename is a default param (nullptr)
//   call with no parameters or nullptr to output
//   to stdout
void
ProblemSpec::output(const char* filename) const
{
  xmlKeepBlanksDefault(0);
  int bytes_written = xmlSaveFormatFileEnc(filename, d_node->doc, "UTF-8", 1);

  if (bytes_written == -1) {
    throw InternalError(std::string("ProblemSpec::output failed for ") +
                          filename,
                        __FILE__,
                        __LINE__);
  }
}

void
ProblemSpec::releaseDocument()
{
  xmlFreeDoc(d_node->doc);
}

ProblemSpecP
ProblemSpec::getRootNode()
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::getRootNode()");
  xmlNode* root_node = xmlDocGetRootElement(d_node->doc);
  return scinew ProblemSpec(
    root_node,
    false); // don't mark as toplevel as this is just a copy
}

const Uintah::TypeDescription*
ProblemSpec::getTypeDescription()
{
  // std::cerr << "ProblemSpec::getTypeDescription() not done\n";
  return 0;
}

// static
ProblemSpecP
ProblemSpec::createDocument(const std::string& name)
{
  MALLOC_TRACE_TAG_SCOPE("ProblemSpec::createDocument()");
  xmlDocPtr doc   = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr node = xmlNewDocRawNode(doc, 0, BAD_CAST name.c_str(), 0);

  xmlDocSetRootElement(doc, node);

  return scinew ProblemSpec(node, true);
}

std::string
ProblemSpec::getFile() const
{
  if (d_node->doc->URL) {
    return (const char*)(d_node->doc->URL);
  } else if (d_node->_private) {
    return (const char*)(d_node->_private);
  } else {
    return "Filename not known.";
  }
}
//______________________________________________________________________
//   Search through the saved labels in the ups:DataArchive section and
//   return true if name is found
bool
ProblemSpec::isLabelSaved(const std::string& name)
{
  ProblemSpecP root  = getRootNode();
  ProblemSpecP DA_ps = root->findBlock("DataArchiver");

  if (!DA_ps) {
    std::string error = "ERROR:  The <DataArchiver> node was not found";
    throw ProblemSetupException(error, __FILE__, __LINE__);
  }

  for (ProblemSpecP var_ps = DA_ps->findBlock("save"); var_ps != 0;
       var_ps              = var_ps->findNextBlock("save")) {

    std::map<std::string, std::string> saveLabel;
    var_ps->getAttributes(saveLabel);
    if (saveLabel["label"] == name) {
      return true;
    }
  }
  return false;
}
