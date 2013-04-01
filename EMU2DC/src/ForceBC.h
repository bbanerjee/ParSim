#ifndef EMU2DC_FORCEBC_H
#define EMU2DC_FORCEBC_H

namespace Emu2DC {
  
  class ForceBC {
  
  public:

   ForceBC();
   ~ForceBC();

   void computeExtForceDensity();

  private:

   // prevent copying
   ForceBC(const ForceBC& dyna);
   ForceBC& operator=(const ForceBC& dyna);

  }; // end class

} // end namespace
#endif

