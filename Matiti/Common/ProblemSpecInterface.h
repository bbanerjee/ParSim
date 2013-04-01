#ifndef MATITI_PROBLEMSPECINTERFACE_H
#define MATITI_PROBLEMSPECINTERFACE_H

#include <Common/SerialPort.h>
#include <Core/ProblemSpec/ProblemSpecP.h>


#include <string>

namespace Matiti {

   class ProblemSpecInterface : public SerialPort {
   public:
      ProblemSpecInterface();
      virtual ~ProblemSpecInterface();

      virtual Uintah::ProblemSpecP readInputFile( const std::string & filename, bool validate = false ) = 0;
      virtual std::string getInputFile() = 0;
      
   private:
      ProblemSpecInterface( const ProblemSpecInterface & );
      ProblemSpecInterface & operator=( const ProblemSpecInterface & );
   };
} // End namespace Matiti

#endif

