#include <Boundary/BoundaryContact.h>
#include <DiscreteElements/DEMParticle.h>

namespace dem {

void 
BoundaryContact::write(std::ostream& os) const
{
  writePosition(os);
  writeNormalForce(os);
  writeTangentForce(os);
  writePenetration(os);
  os << std::endl;
}

void 
BoundaryContact::writePosition(std::ostream& os) const
{
  os << std::setw(OWID) << bc_point.x() 
     << std::setw(OWID) << bc_point.y()
     << std::setw(OWID) << bc_point.z() ;
}

void 
BoundaryContact::writeNormalForce(std::ostream& os) const
{
  os << std::setw(OWID) << bc_normalForce.x()
     << std::setw(OWID) << bc_normalForce.y() 
     << std::setw(OWID) << bc_normalForce.z();
}

void 
BoundaryContact::writeTangentForce(std::ostream& os) const
{
  os << std::setw(OWID) << bc_tangentForce.x()
     << std::setw(OWID) << bc_tangentForce.y() 
     << std::setw(OWID) << bc_tangentForce.z();
}

void 
BoundaryContact::writePenetration(std::ostream& os) const
{
  os << std::setw(OWID) << bc_penetration;
}

void 
BoundaryContact::write(zen::XmlOut& xml) const
{
  writePosition(xml);
  writeNormalForce(xml);
  writeTangentForce(xml);
  writePenetration(xml);
}

void 
BoundaryContact::writePosition(zen::XmlOut& xml) const
{
  std::ostringstream stream;
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.precision(dem::OPREC);
  stream << "["
        << bc_point.x() << ", " << bc_point.y() << ", " << bc_point.z()
        << "]";
  xml["position"](stream.str());
}

void 
BoundaryContact::writeNormalForce(zen::XmlOut& xml) const
{
  std::ostringstream stream;
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.precision(dem::OPREC);
  stream << "["
        << bc_normalForce.x() << ", " 
        << bc_normalForce.y() << ", " 
        << bc_normalForce.z()
        << "]";
  xml["normal_force"](stream.str());
}

void 
BoundaryContact::writeTangentForce(zen::XmlOut& xml) const
{
  std::ostringstream stream;
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.precision(dem::OPREC);
  stream << "["
        << bc_tangentForce.x() << ", " 
        << bc_tangentForce.y() << ", " 
        << bc_tangentForce.z()
        << "]";
  xml["tangent_force"](stream.str());
}

void 
BoundaryContact::writePenetration(zen::XmlOut& xml) const
{
  xml["penetration"](bc_penetration);
}

} // end namespace dem