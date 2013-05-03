#ifndef EMU2DC_LINE_H
#define EMU2DC_LINE_H

namespace Emu2DC {
    
  class Line {
    
    public:
      Line();
      virtual ~Line();

      bool intersects(const Line& line);

    protected:

      double x1, y1, x2, y2;
    
  }; // end class
} // end namespace

#endif

