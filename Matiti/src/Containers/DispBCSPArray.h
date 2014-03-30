#ifndef MATITI_DISP_BC_SP_ARRAY_H
#define MATITI_DISP_BC_SP_ARRAY_H

#include <Pointers/DispBCSP.h>
#include <vector>

namespace Matiti {
  
  typedef std::vector<DispBCSP> DispBCSPArray;
  typedef std::vector<DispBCSP>::iterator DispBCSPIterator;
  typedef std::vector<DispBCSP>::const_iterator constDispBCSPIterator;
}

#endif
