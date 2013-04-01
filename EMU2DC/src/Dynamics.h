#ifndef EMU2DC_DYNAMICS_H
#define EMU2DC_DYNAMICS_H

//********************************************************************
// ModuleFile : dynamic_integration
// Purpose    : This module includes subroutines to calculate the 
//              internal and external force density, which are used to 
//              calculate displacements
//********************************************************************

namespace Emu2DC {
  
  class Dynamics {
  
  public:

   Dynamics();
   ~Dynamics();

   void do_dynamics();
   void write_output(const int& output_file_count);

  private:

   // prevent copying
   Dynamics(const Dynamics& dyna);
   Dynamics& operator=(const Dynamics& dyna);

  }; // end class

} // end namespace
#endif

