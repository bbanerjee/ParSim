/*
 * MPMConstitutiveModel.cc
 *
 *  Created on: 21/10/2013
 *      Author: banerjee
 */

#include "MPMConstitutiveModel.h"

using namespace BrMPM;

MPMConstitutiveModel::MPMConstitutiveModel() {
	// TODO Auto-generated constructor stub

}

MPMConstitutiveModel::~MPMConstitutiveModel() {
	// TODO Auto-generated destructor stub
}

void BrMPM::MPMConstitutiveModel::getStress(const Matrix3D& defGrad,
                                            Matrix3D& stress,
                                            double& jacobian)
{
  // TODO  Hooman to add in details.  Use a factory.
}

