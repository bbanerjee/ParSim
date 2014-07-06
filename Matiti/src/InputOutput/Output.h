#ifndef __MATITI_OUTPUT_H__
#define __MATITI_OUTPUT_H__

#include <Core/Time.h>
#include <Core/Domain.h>
#include <Containers/BodySPArray.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <string>
#include <iostream>

namespace Matiti {

  class Output 
  {
  public:
    friend std::ostream& operator<<(std::ostream& out, const Matiti::Output& output);

  public:
    Output();
    Output(const std::string& fileName,
           int iterInterval);
    Output(const Uintah::ProblemSpecP& ps);
    virtual ~Output();

    void clone(const Output& output);

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

