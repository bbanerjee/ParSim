/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2014-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

#ifndef __PROBLEM_SPEC_READER_H__
#define __PROBLEM_SPEC_READER_H__

#include <CCA/Ports/ProblemSpecInterface.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <vector>
#include <memory>

namespace Uintah {

struct Tag;

class UPSFilenameManager
{
  std::vector<std::unique_ptr<std::string>> filenames_;

public:
  std::string*
  store(const std::string& filename)
  {
    filenames_.emplace_back(std::make_unique<std::string>(filename));
    return filenames_.back().get();
  }

  std::string
  getMainXMLFilename() const
  {
    return *(filenames_[0]);
  }

  // Find a filename in the stored collection
  const std::string&
  findFileName(const std::string& filename) const
  {
    for (const auto& stored : filenames_) {
      if (*stored == filename) {
        return *stored;
      }
    }
    static const std::string none = "None";
    return none;
  }

  // Non-const version that returns a pointer (for _private assignment)
  std::string*
  findFileNamePtr(const std::string& filename)
  {
    for (const auto& stored : filenames_) {
      if (*stored == filename) {
        return stored.get();
      }
    }
    return nullptr;
  }

  // Store if not found, otherwise return existing
  std::string*
  storeOrFind(const std::string& filename)
  {
    if (auto* existing = findFileNamePtr(filename)) {
      return existing;
    }
    return store(filename);
  }

  // Get the number of stored filenames
  size_t
  size() const
  {
    return filenames_.size();
  }

  // Check if empty
  bool
  empty() const
  {
    return filenames_.empty();
  }

  // Clear all stored filenames
  void
  clear()
  {
    filenames_.clear();
  }

  // Automatic cleanup when FilenameManager is destroyed
  ~UPSFilenameManager() = default;

  // Automatic cleanup when UPSFilenameManager is destroyed
};

class ProblemSpecReader : public ProblemSpecInterface
{
public:
  ProblemSpecReader() = default;
  ~ProblemSpecReader();

  // Be sure to call releaseDocument on this ProblemSpecP.  Most users should
  // not use the 'insertAfterThisNode' parameter.
  ProblemSpecP readInputFile(const std::string& filename,
                             bool validate = false) override;

  // Returns the main xml file name.
  std::string getInputFile() override { return d_upsFilename.getMainXMLFilename(); }

  // Set the xml data and file name
  // ** WARNING ** Only for testing purposes
  void setData(ProblemSpecP data) override;
  void setFilename(const std::string& filename) override;

private:

  ProblemSpecReader(const ProblemSpecReader&) = delete;
  ProblemSpecReader& operator=(const ProblemSpecReader&) = delete;

  ////////////////////////////////////////////////////////////////////////////////
  // Variables:

  // d_upsFilename[0] is the main file... but each subsequent string
  // is the name of an <include>d file.
  UPSFilenameManager d_upsFilename;

  ProblemSpecP d_xmlData;

  ////////////////////////////////////////////////////////////////////////////////
  // Functions:

  void parseValidationFile();
  void validateProblemSpec(ProblemSpecP& prob_spec);
  std::string* findFileNamePtr(const std::string& filename);

  // Replaces <include> tags with xml file tree.
  void resolveIncludes(xmlNode* child, xmlNode* parent, int depth = 0);
};

} // End namespace Uintah

#endif // __PROBLEM_SPEC_READER_H__
