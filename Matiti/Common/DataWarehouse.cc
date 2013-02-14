#include <Common/DataWarehouse.h>
#include <Common/Scheduler.h>

#include <iostream>

using namespace Matiti;

DataWarehouse::DataWarehouse(Scheduler* scheduler,
			     int generation)
  : d_scheduler(scheduler), d_generation(generation)
{
}

DataWarehouse::~DataWarehouse()
{
}
