#ifndef EMU2DC_DISPBC_H
#define EMU2DC_DISPBC_H

namespace Emu2DC {
  
  class DispBC {
  
  public:

   DispBC();
   ~DispBC();

   void apply_boundary_cond();

  private:

   // prevent copying
   DispBC(const DispBC& bc);
   DispBC& operator=(const DispBC& bc);

  }; // end class

} // end namespace

#endif

