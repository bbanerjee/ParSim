#ifndef EMU2DC_INITIALCONDITIONS_H
#define EMU2DC_INITIALCONDITIONS_H

namespace Emu2DC {
  
  class InitialConditions {
  
  public:

   InitialConditions();
   ~InitialConditions();

   void apply_initial_conditions();

  private:

   // prevent copying
   InitialConditions(const InitialConditions& bc);
   InitialConditions& operator=(const InitialConditions& bc);

  }; // end class

} // end namespace
#endif

