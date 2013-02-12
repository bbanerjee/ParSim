#include <Common/SerialComponent.h>
#include <algorithm>

using namespace Matiti;
using std::map;
using std::string;

SerialComponent::SerialComponent()
{
}

SerialComponent::~SerialComponent()
{
  for(map<string, PortRecord*>::iterator iter = portmap.begin(); 
      iter != portmap.end(); iter++) 
    delete iter->second;

}

void
SerialComponent::attachPort(const string& name,
                            SerialPort* port)
{
    map<string, PortRecord*>::iterator iter = portmap.find(name);
    if(iter == portmap.end()){
      portmap[name]=scinew PortRecord(port);
    } else {
      iter->second->connections.push_back(port);
    }
}

SerialComponent::PortRecord::PortRecord(SerialPort* port)
{
    connections.push_back(port);
}

SerialPort* SerialComponent::getPort(const std::string& name)
{
    map<string, PortRecord*>::iterator iter = portmap.find(name);
    if(iter == portmap.end())
      return 0;
    else if(iter->second->connections.size()> 1)
      return iter->second->connections.back();
    else
      return iter->second->connections[0];
}

SerialPort* SerialComponent::getPort(const std::string& name,
                                     unsigned int i)
{
    map<string, PortRecord*>::iterator iter = portmap.find(name);
    if(iter == portmap.end())
      return 0;
    else if(iter->second->connections.size()> 1)
      return iter->second->connections[i];
    else
      return iter->second->connections[0];
}

void SerialComponent::releasePort(const std::string&)
{
}

unsigned int SerialComponent::numConnections(const std::string& name)
{
  map<string, PortRecord*>::iterator iter = portmap.find(name);
  if(iter == portmap.end())
    return 0;
  else 
    return iter->second->connections.size();
}
