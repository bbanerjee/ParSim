#ifndef MATITI_OUTPUTCONTEXT_H
#define MATITI_OUTPUTCONTEXT_H

#include <Core/ProblemSpec/ProblemSpec.h>

namespace Matiti {
   class OutputContext {
   public:
      OutputContext(int fd, const char* filename, long cur, Uintah::ProblemSpecP varnode)
	: fd(fd), filename(filename), cur(cur), varnode(varnode)
      {
      }
      ~OutputContext() {}

      int fd;
      const char* filename;
      long cur;
      Uintah::ProblemSpecP varnode;
   private:
      OutputContext(const OutputContext&);
      OutputContext& operator=(const OutputContext&);
      
   };
} // End namespace Matiti

#endif
