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
#include <Matrix3D.h>

#include <vector>
#include <string>

namespace BrMPM {

  template <class T> class MPMData {

    public:

      MPMData() {}
      virtual ~MPMData() {}

      typename std::vector<T>::iterator begin() {return d_data.begin();}
      typename std::vector<T>::iterator end() {return d_data.end();}
      typename std::vector<T>::const_iterator begin() const {return d_data.begin();}
      typename std::vector<T>::const_iterator end() const {return d_data.end();}

      unsigned int size() {return d_data.size();}
      unsigned int size() const {return d_data.size();}

      T& operator[](int index) {return d_data.at(index);}
      const T& operator[](int index) const {return d_data.at(index);}

    private:

      int d_material_index;
      std::string d_variable_name;
      std::vector<T> d_data;

  }; // end class

  // Define instances
  // **WARNING**  These are per particle variables.  A ParticleData object contains
  // values for all particles
  typedef MPMData<int>                     IntegerParticleData;
  typedef MPMData<double>                  DoubleParticleData;
  typedef MPMData<std::string>             StringParticleData;
  typedef MPMData<Point3D>                 Point3DParticleData;
  typedef MPMData<Vector3D>                Vector3DParticleData;
  typedef MPMData<Matrix3D>                Matrix3DParticleData;
  typedef MPMData<std::vector<int> >       VectorIntParticleData;
  typedef MPMData<std::vector<double> >    VectorDoubleParticleData;
  typedef MPMData<std::vector<Vector3D> >  VectorVector3DParticleData;
  typedef MPMData<std::vector<Matrix3D> >  VectorMatrix3DParticleData;

  typedef MPMData<double>                  DoubleNodeData;
  typedef MPMData<Vector3D>                Vector3DNodeData;

} /* namespace BrMPM */

#endif /* MPMDATA_H_ */
