/*
 * MPMConstitutiveModel.h
 *
 *  Created on: 21/10/2013
 *      Author: banerjee
 */

#ifndef MPMCONSTITUTIVEMODEL_H_
#define MPMCONSTITUTIVEMODEL_H_

#include <GeometryMath/Matrix3D.h>

namespace BrMPM {

class MPMConstitutiveModel {
public:
	MPMConstitutiveModel();
	virtual ~MPMConstitutiveModel();

	void getStress(const Matrix3D& defGrad,
	               Matrix3D& stress,
	               double& jacobian);
};

} /* namespace BrMPM */
#endif /* MPMCONSTITUTIVEMODEL_H_ */
