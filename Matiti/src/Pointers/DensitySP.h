#ifndef MATITI_DENSITY_SP_H
#define MATITI_DENSITY_SP_H

#include <memory>

namespace Matiti {

 // Forward declaration, Make sure <Density.h> is included before using DensitySP.
 // Using stdlib shared_ptr

 class Density;
 typedef std::shared_ptr<Density> DensitySP;
 }

#endif 
