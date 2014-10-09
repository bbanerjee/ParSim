/*
 * MPMContactFactory.h
 *
 *  Created on: 23/10/2013
 *      Author: banerjee
 */

#ifndef MPMCONTACTFACTORY_H_
#define MPMCONTACTFACTORY_H_

#include <Contact/MPMContactP.h>
#include <MPMPatchP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <vector>

namespace BrMPM
{

class MPMContactFactory
{
public:
  MPMContactFactory();
  virtual ~MPMContactFactory();
  static MPMContactP create(const Uintah::ProblemSpecP& ps,
                            std::vector<int>& dwis,
                            MPMPatchP& patch);
};

} /* namespace BrMPM */
#endif /* MPMCONTACTFACTORY_H_ */
