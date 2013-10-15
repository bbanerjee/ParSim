/*
 * MPMData.h
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#ifndef MPMDATA_H_
#define MPMDATA_H_

#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>
#include <MPMMatrix.h>

#include <vector>
#include <string>

namespace BrMPM {

  template <class T> class MPMData {

    public:

      MPMData() {}
      virtual ~MPMData() {}

      std::vector<T>::iterator begin() {return d_data.begin();}
      std::vector<T>::iterator end() {return d_data.end();}

      unsigned int size() {return d_data.size();}

    private:

      std::vector<T> d_data;

  }; // end class

  // Define instances
  // **WARNING**  These are per particle variables.  A ParticleData object contains
  // values for all particles
  typedef MPMData<int>                   IntegerParticleData;
  typedef MPMData<double>                DoubleParticleData;
  typedef MPMData<std::string>           StringParticleData;
  typedef MPMData<Point3D>               PointParticleData;
  typedef MPMData<Vector3D>              VectorParticleData;
  typedef MPMData<Matrix3D>              MatrixParticleData;
  typedef MPMData<std::vector<int> >     VectorIntParticleData;
  typedef MPMData<std::vector<double> >  VectorDoubleParticleData;
  typedef MPMData<double>                DoubleNodeData;

} /* namespace BrMPM */

#endif /* MPMDATA_H_ */
