/*
 * ExtensionVelocityMarcher.h
 *
 *  Created on: 5/12/2013
 *      Author: banerjee
 */

#ifndef EXTENSIONVELOCITYMARCHER_H_
#define EXTENSIONVELOCITYMARCHER_H_

#include <Contact/FastMarching/DistanceMarcher.h>

namespace BrMPM
{

  class ExtensionVelocityMarcher: public DistanceMarcher
  {
  public:
    ExtensionVelocityMarcher(double* phi, double* dx, long* flag,
                      double* distance, int ndim, const std::vector<int>& shape,
                      bool self_test, int order,
                      double* speed, double* f_ext, long* ext_mask);
    virtual ~ExtensionVelocityMarcher();

  protected:

    virtual void initializeFrozen();
    virtual void finalizePoint(int ii, double phi_i);
    virtual void cleanUp();

  private:

    double* d_speed;
    double* d_f_ext;
    long* d_ext_mask;
  };

} /* namespace BrMPM */
#endif /* EXTENSIONVELOCITYMARCHER_H_ */
