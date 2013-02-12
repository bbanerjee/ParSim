#ifndef MATITI_SERIALPORT_H
#define MATITI_SERIALPORT_H

#include <Common/RefCounted.h>


namespace Matiti {


class SerialPort : public RefCounted {
public:
           SerialPort();
  virtual ~SerialPort();
};

} // End namespace Matiti

#endif
