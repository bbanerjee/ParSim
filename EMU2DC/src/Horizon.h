#ifndef EMU2DC_HORIZON_H
#define EMU2DC_HORIZON_H

namespace Emu2DC {

  class Horizon {

  public:

    Horizon();
    ~Horizon();

    //****************************************************************************************************
    // This subroutine will calculate the horizon)size for each node of any uniform and non-uniform grid. 
    // This value is the largest edge connceted to this node times an scalar defined by the user
    //!*****************************************************************************************************
    void calculateHorizon();

  private:

    // Prevent copying
    Horizon(const Horizon& horizon);
    Horizon& operator=(const Horizon& horizon);

  }; // end class Horizon
}  // End namespace

#endif
