#ifndef MATITI_WOOD_SP_H
#define MATITI_WOOD_SP_H

#include <memory>

namespace Matiti {

 // Forward declaration, Make sure <Wood.h> is included before using DensitySP.
 // Using stdlib shared_ptr

 class Wood;
 typedef std::shared_ptr<Wood> WoodSP;
 }

#endif 
