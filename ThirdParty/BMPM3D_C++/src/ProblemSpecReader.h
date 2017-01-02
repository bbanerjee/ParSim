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

#ifndef __MATITI_PROBLEM_SPEC_READER_H__ 
#define __MATITI_PROBLEM_SPEC_READER_H__

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <string>
#include <vector>

namespace BrMPM {
      
  class ProblemSpecReader {

  public:
    ProblemSpecReader();
    ~ProblemSpecReader();

    // Be sure to call releaseDocument on this ProblemSpecP.  
    virtual Uintah::ProblemSpecP readInputFile( const std::string & filename);

    // Returns the main xml file name.
    virtual std::string getInputFile() { return *d_upsFilename[0]; }

  private:

    ProblemSpecReader(const ProblemSpecReader&);
    ProblemSpecReader& operator=(const ProblemSpecReader&);
    
    ////////////////////////////////////////////////////////////////////////////////
    // Variables:

    // d_upsFilename[0] is the main file... but each subsequent string
    // is the name of an <include>d file.
    std::vector< std::string * > d_upsFilename;

    Uintah::ProblemSpecP d_xmlData;
   };

} // End namespace 

#endif // __PROBLEM_SPEC_READER_H__
