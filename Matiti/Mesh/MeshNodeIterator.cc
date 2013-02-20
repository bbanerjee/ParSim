#include <Mesh/MeshNodeIterator.h>


namespace Matiti
{
  std::ostream& operator<<(std::ostream& out, const Matiti::MeshNodeIterator& b)
  {
    out << "[MeshNodeIterator at " << *b << " of " << b.end() << ']';

    return out;
  }
}
