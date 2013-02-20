#include <Mesh/MeshElementIterator.h>
#include <iostream>

using namespace Matiti;


namespace Matiti
{
  ostream& operator<<(ostream& out, const MeshElementIterator& c)
  {
    out << "[MeshElementIterator at " << *c << " of " << c.end() << ']';
    return out;
  }
}

