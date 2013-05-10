#ifndef __EMU2DC_OUTPUT_H__
#define __EMU2DC_OUTPUT_H__

#include <Time.h>
#include <BodySPArray.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <iostream>

namespace Emu2DC {

  class Output 
  {
  public:
    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::Output& output);

  public:
    Output();
    Output(const Uintah::ProblemSpecP& ps);
    ~Output();

    void initialize(const Uintah::ProblemSpecP& ps);
    void write(const Time& time, const BodySPArray& bodyList);

    inline void outputFolder(const std::string& folder) {d_output_folder_name = folder;}
    inline std::string outputFolder() const {return d_output_folder_name;}
    inline std::string outputFile() const {return d_output_file_name;}
    inline int outputIteratonInterval() const {return d_output_iter_interval;}

    int outputFileCount() const {return d_output_file_count;}

  private:

    void incrementOutputFileCount() {d_output_file_count++;}

    //  Output file folder and name 
    std::string d_output_folder_name;
    std::string d_output_file_name;
    int d_output_iter_interval;

    int d_output_file_count;
  }; // end class

} // end namespace

#endif

