/*
 * MPMContact.cc
 *
 *  Created on: 14/10/2013
 *      Author: banerjee
 */

#include "MPMContact.h"

using namespace BrMPM;

MPMContact::MPMContact(std::vector<int>& dwis,
                       MPMPatchP& patch)
  : d_dwis(dwis), d_patch(patch), d_mtol(1.0e-10)
{
}

MPMContact::~MPMContact()
{
}

void
MPMContact::findIntersection(MPMDatawarehouseP& dw)
{
}

void
MPMContact::findIntersectionSimple(MPMDatawarehouseP& dw)
{
}
