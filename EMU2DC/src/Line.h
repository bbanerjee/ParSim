#ifndef EMU2DC_LINE
#define EMU2DC_LINE

namespace Emu2DC {
    
  class Line {
    
    public:
      Line();
      ~Line();
      bool intersects(const Line& line);

    protected:

      double x1, y1, x2, y2;
    
  };
} // end namespace

#endif

