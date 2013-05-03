#ifndef EMU2DC_HORIZON_COMPUTER_H
#define EMU2DC_HORIZON_COMPUTER_H

#include <BodySP.h>
#include <SimulationState.h>

//************************************** 
/** 
  * @file
  * @author  Biswajit Banerjee (2013)
  * @version 1.0
  *
  * @section LICENSE
  *
  * Copyright (c) 2013 Callaghan Innovation
  *
  * @section DESCRIPTION
  *
  * This functor will calculate the horizon size for each node of any uniform and non-uniform grid. 
  * This value is the largest edge connceted to this node times an scalar defined by the user
  */

namespace Emu2DC {

  class HorizonComputer {

  public:

    HorizonComputer() {}

   /**
    * Compute the horizons of all the nodes in the body
    *
    * Usage:
    *
    *  HorizonComputer compute_horizon;
    *  compute_horizon(body, state);
    */
    void operator()(const BodySP& body, SimulationState& state);

  }; // end class HorizonComputer
}  // End namespace

#endif
