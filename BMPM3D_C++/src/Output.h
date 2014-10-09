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

#ifndef __MATITI_OUTPUT_H__
#define __MATITI_OUTPUT_H__

#include <Time.h>
#include <Domain.h>
#include <BodySPArray.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <iostream>

namespace BrMPM {

  class Output 
  {
  public:
    friend std::ostream& operator<<(std::ostream& out, const BrMPM::Output& output);

  public:
    Output();
    Output(const Uintah::ProblemSpecP& ps);
    virtual ~Output();

    void initialize(const Uintah::ProblemSpecP& ps);
    virtual void write(const Time& time, const Domain& domain, const BodySPArray& bodyList);

    inline void outputFolder(const std::string& folder) {d_output_folder_name = folder;}
    inline std::string outputFolder() const {return d_output_folder_name;}
    inline std::string outputFile() const {return d_output_file_name;}
    inline int outputIteratonInterval() const {return d_output_iter_interval;}

    int outputFileCount() const {return d_output_file_count;}

  protected:

    void incrementOutputFileCount() {d_output_file_count++;}

  private:

    //  Output file folder and name 
    std::string d_output_folder_name;
    std::string d_output_file_name;
    int d_output_iter_interval;

    int d_output_file_count;
  }; // end class

} // end namespace

#endif

