/*
 * MPMContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include <Contact/MPMContact.h>
#include <MPMDatawarehouse.h>
#include <MPMPatch.h>

using namespace BrMPM;

// Initialize default contact algorithm
MPMContact::MPMContact()
  :d_mtol(1.0e-10)
{
}

// Initialize default contact algorithm
MPMContact::MPMContact(std::vector<int>& dwis,
                       MPMPatchP& patch)
  : d_dwis(dwis), d_patch(patch), d_mtol(1.0e-10)
{
}

MPMContact::~MPMContact()
{
}

// Find whether two bodies intersect using a fast marching algorithm
void
MPMContact::findIntersection(MPMDatawarehouseP& dw)
{

}

void
MPMContact::findIntersectionSimple(MPMDatawarehouseP& dw)
{
}

void
MPMContact::exchMomentumInterpolated(MPMDatawarehouseP& dw)
{
}

void
MPMContact::exchForceInterpolated(MPMDatawarehouseP& dw)
{
}

void
MPMContact::exchMomentumIntegrated(MPMDatawarehouseP& dw)
{
}
