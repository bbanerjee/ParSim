#ifndef MATITI_SERIALCOMPONENT_H
#define MATITI_SERIALCOMPONENT_H

#include <string>
#include <map>
#include <vector>

namespace Matiti {

  using std::string;
  class SerialPort;

  class SerialComponent {
    struct PortRecord {
      PortRecord(SerialPort* conn);
      std::vector<SerialPort*> connections;
    };
    std::map<string, PortRecord*> portmap;

    public:
      SerialComponent();
      virtual ~SerialComponent();
      
      void attachPort(const string& name, SerialPort* port);
      
      SerialPort* getPort(const std::string& name);
      SerialPort* getPort(const std::string& name, unsigned int i);
      void releasePort(const std::string& name);
      unsigned int numConnections(const std::string& name);
   };
} // End namespace Matiti
   
#endif
