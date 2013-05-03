#ifndef __EMU2DC_PROBLEM_SPEC_READER_H__ 
#define __EMU2DC_PROBLEM_SPEC_READER_H__

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <string>
#include <vector>

namespace Emu2DC {
      
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
