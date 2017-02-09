#ifndef BOUNDARY_CONTAINERS_H
#define BOUNDARY_CONTAINERS_H

#include <memory>
#include <vector>

namespace dem {

class Boundary;
class BoundaryTangent;
class BoundaryContact;

using BoundaryUP = std::unique_ptr<Boundary>;
using BoundaryUPArray = std::vector<BoundaryUP>;

using BoundaryP = std::shared_ptr<Boundary>;
using BoundaryPArray = std::vector<BoundaryP>;

using BoundaryTangentArray = std::vector<BoundaryTangent>;
using BoundaryTangentArrayMap = std::map<std::size_t, BoundaryTangentArray>;

using BoundaryContactArray = std::vector<BoundaryContact>;
}

#endif
