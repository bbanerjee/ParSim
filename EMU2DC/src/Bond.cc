#include <Bond.h>

using namespace Emu2DC;

Bond::Bond():d_node1(0),d_node2(0),d_broken(false)
{
}

Bond::Bond(Node* node1, Node* node2):d_node1(node1),d_node2(node2),d_broken(false)
{
}
Bond::~Bond()
{
}

void 
Bond::first(const Node* node)
{
  d_node1 = node;
}

void 
Bond::second(const Node* node)
{
  d_node2 = node;
}

void 
Bond::isBroken(const bool broken)
{
  d_broken = broken;
}

Node* 
Bond::first() const
{
  return d_node1;
}

Node* 
Bond::second() const
{
  return d_node2;
}

bool 
Bond::isBroken() const
{
  return d_broken;
}
    
bool 
Bond::operator==(const Bond& bond) const
{
  return ((d_node1 == bond.d_node1 && d_node2 == bond.d_node2) ||
          (d_node1 == bond.d_node2 && d_node2 == bond.d_node1));
}

