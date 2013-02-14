#include <Common/ComputeSet.h>
#include <iostream>


using namespace Matiti;

namespace Matiti {

  ostream& operator<<(ostream& out, const Matiti::MaterialSubset& mss)
  {
    out << "{";
    for(int j=0;j<mss.size();j++){
      if(j != 0)
        out << ",";
      out << mss.get(j);
    }
    out << "}";
    return out;
  }

  ostream& operator<<(ostream& out, const Matiti::MaterialSet& ms)
  {
    if(&ms == 0)
      out << "(null Materials)";
    else {
      out << "Matls: {";
      for(int i=0;i< ms.size();i++){
        const MaterialSubset* mss = ms.getSubset(i);
        if(i != 0)
          out << ", ";
        out << *mss;
      }
      out << "}";
    }
    return out;
  }

} // end namespace Matiti


  
