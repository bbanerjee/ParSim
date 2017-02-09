#include <Boundary/Boundary.h>
#include <DiscreteElements/Particle.h>
// use both pointer to and variable of class Particle

namespace dem {

Boundary::Boundary(std::size_t tp, std::ifstream& ifs)
{
  type = tp;
  ifs >> extraNum;
  ifs >> id;
}

void
Boundary::print(std::ostream& os)
{
  os << std::endl
     << std::setw(OWID) << type << std::setw(OWID) << extraNum << std::endl
     << std::setw(OWID) << id;
}

void
Boundary::printContactInfo(std::ostream& os)
{
  os << std::setw(OWID) << id << std::endl
     << std::setw(OWID) << contactInfo.size() << std::endl
     << std::setw(OWID) << "pos_x" << std::setw(OWID) << "pos_y"
     << std::setw(OWID) << "pos_z" << std::setw(OWID) << "normal_x"
     << std::setw(OWID) << "normal_y" << std::setw(OWID) << "normal_z"
     << std::setw(OWID) << "tangt_x" << std::setw(OWID) << "tangt_y"
     << std::setw(OWID) << "tangt_z" << std::setw(OWID) << "pentr" << std::endl;

  for (BoundaryContactArray::iterator it = contactInfo.begin();
       it != contactInfo.end(); ++it)
    it->print(os);
}

void
Boundary::clearStatForce()
{
  contactNum = 0;
  normal = 0;
  tangt = 0;
  penetr = 0;
}

void
Boundary::updateStatForce()
{
  clearStatForce();
  contactNum = contactInfo.size();
  for (BoundaryContactArray::iterator it = contactInfo.begin();
       it != contactInfo.end(); ++it) {
    normal += it->normal;
    tangt += it->tangt;
    penetr += it->penetr;
  }
  if (contactNum != 0)
    penetr /= contactNum;
}

void
Boundary::clearContactInfo()
{
  possParticle.clear();
  contactInfo.clear();
}

} // namespace dem ends
