#ifndef BOUNDARY_CONTACT_H
#define BOUNDARY_CONTACT_H

#include <Core/Const/Constants.h>
#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>
#include <boost/serialization/base_object.hpp>
#include <InputOutput/zenxml/xml.h>
#include <iostream>
#include <utility>

namespace dem {

class DEMParticle; // forward declaration, only use pointer to class DEMParticle

/////////////////////////////////////
class BoundaryContact
{
public:
  DEMParticle* bc_particle;
  Vec bc_point;
  Vec bc_normalForce;
  Vec bc_tangentForce;
  REAL bc_penetration;

public:
  BoundaryContact()
    : bc_particle(nullptr)
    , bc_point(0)
    , bc_normalForce(0)
    , bc_tangentForce(0)
    , bc_penetration(0)
  {
  }

  BoundaryContact(DEMParticle* p, Vec pt, Vec nm, Vec tg, REAL pntr)
    : bc_particle(p)
    , bc_point(std::move(pt))
    , bc_normalForce(std::move(nm))
    , bc_tangentForce(std::move(tg))
    , bc_penetration(pntr)
  {
  }

  void write(std::ostream& os) const;
  void write(zen::XmlOut& xml) const;

private:

  inline void writePosition(std::ostream& os) const;
  inline void writeNormalForce(std::ostream& os) const;
  inline void writeTangentForce(std::ostream& os) const;
  inline void writePenetration(std::ostream& os) const;

  void writePosition(zen::XmlOut& xml) const;
  void writeNormalForce(zen::XmlOut& xml) const;
  void writeTangentForce(zen::XmlOut& xml) const;
  void writePenetration(zen::XmlOut& xml) const;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& bc_particle;
    ar& bc_point;
    ar& bc_normalForce;
    ar& bc_tangentForce;
    ar& bc_penetration;
  }
};

} // end namespace dem

#endif
