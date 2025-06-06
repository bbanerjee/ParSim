/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <CCA/Components/LoadBalancers/ProfileDriver.h>
#include <Core/Util/DebugStream.h>
#include <iomanip>
#include <sstream>

using namespace Uintah;

static DebugStream stats("ProfileStats", false);
static DebugStream stats2("ProfileStats2", false);

//______________________________________________________________________
//
void
ProfileDriver::setMinPatchSize(const std::vector<IntVector>& min_patch_size)
{
  d_minPatchSize = min_patch_size;
  costs.resize(d_minPatchSize.size());
  d_minPatchSizeVolume.resize(d_minPatchSize.size());
  for (int l = 0; l < (int)d_minPatchSize.size(); l++) {
    d_minPatchSizeVolume[l] =
      d_minPatchSize[l][0] * d_minPatchSize[l][1] * d_minPatchSize[l][2];
  }
}
//______________________________________________________________________
//
void
ProfileDriver::addContribution(const PatchSubset* patches, double cost)
{
  if (patches) {
    // compute number of datapoints
    int num_points = 0;
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      int l              = patch->getLevel()->getIndex();
      num_points += patch->getNumCells() / d_minPatchSizeVolume[l];
    }
    double average_cost = cost / num_points;

    // loop through patches
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      int l              = patch->getLevel()->getIndex();

      // coarsen region by minimum patch size
      IntVector low  = patch->getCellLowIndex() / d_minPatchSize[l];
      IntVector high = patch->getCellHighIndex() / d_minPatchSize[l];

      // loop through datapoints
      for (CellIterator iter(low, high); !iter.done(); iter++) {
        // add cost to current contribution
        costs[l][*iter].current += average_cost;
        // if(d_myworld->myRank()==0)
        //    std::cout << "level:" << l << " adding: " << average_cost << " to:
        //    " << *iter << " new total:" << costs[l][*iter].current << std::endl;
      }
    }
  } else {
    // if(d_myworld->myRank()==0)
    //    std::cout << "no patches\n";
  }
}

//______________________________________________________________________
//
void
ProfileDriver::outputError(const GridP currentGrid)
{
  static int iter = 0;
  iter++;
  std::vector<double> proc_costsm(d_myworld->nRanks(), 0);
  std::vector<double> proc_costsp(d_myworld->nRanks(), 0);
  std::vector<std::vector<double>> predicted_sum(currentGrid->numLevels()),
    measured_sum(currentGrid->numLevels());

  double smape = 0, max_error = 0, smpe = 0;
  int num_patches = 0;

  //__________________________________
  // each level
  for (int l = 0; l < currentGrid->numLevels(); l++) {
    LevelP level = currentGrid->getLevel(l);
    num_patches += level->numPatches();
    std::vector<Region> regions(level->numPatches());

    predicted_sum[l].assign(level->numPatches(), 0);
    measured_sum[l].assign(level->numPatches(), 0);

    for (int p = 0; p < level->numPatches(); p++) {
      const Patch* patch = level->getPatch(p);
      regions[p] = Region(patch->getCellLowIndex(), patch->getCellHighIndex());
    }

    std::vector<double> predicted(regions.size(), 0),
      measured(regions.size(), 0);

    for (int r = 0; r < (int)regions.size(); r++) {
      // coarsen region by minimum patch size
      IntVector low  = regions[r].getLow() / d_minPatchSize[l];
      IntVector high = regions[r].getHigh() / d_minPatchSize[l];

      // loop through datapoints
      for (CellIterator iter(low, high); !iter.done(); iter++) {
        // search for point in the map
        std::map<IntVector, Contribution>::iterator it = costs[l].find(*iter);

        // if in the map
        if (it != costs[l].end()) {
          // add predicted and measured costs to respective arrays
          predicted[r] += it->second.weight;
          measured[r] += it->second.current;
        }
      }
    }

    //__________________________________
    // allreduce sum weights
    if (d_myworld->nRanks() > 1) {
      Uintah::MPI::Allreduce(&predicted[0],
                             &predicted_sum[l][0],
                             predicted.size(),
                             MPI_DOUBLE,
                             MPI_SUM,
                             d_myworld->getComm());
      Uintah::MPI::Allreduce(&measured[0],
                             &measured_sum[l][0],
                             measured.size(),
                             MPI_DOUBLE,
                             MPI_SUM,
                             d_myworld->getComm());
    }

    //__________________________________
    //
    for (int p = 0; p < level->numPatches(); p++) {
      const Patch* patch = level->getPatch(p);

      int proc = d_lb->getPatchwiseProcessorAssignment(patch);
      proc_costsm[proc] += measured_sum[l][p];
      proc_costsp[proc] += predicted_sum[l][p];
      double error = (measured_sum[l][p] - predicted_sum[l][p]) /
                     (measured_sum[l][p] + predicted_sum[l][p]);

      if (fabs(error) > max_error) {
        max_error = fabs(error);
      }
      smape += fabs(error);
      smpe += error;

      if (d_myworld->myRank() == 0 && stats2.active()) {
        IntVector low(patch->getCellLowIndex()),
          high(patch->getCellHighIndex());
        stats2 << iter << " " << patch->getID() << " "
               << (measured_sum[l][p] - predicted_sum[l][p]) /
                    (measured_sum[l][p] + predicted_sum[l][p])
               << " " << measured_sum[l][p] << " " << predicted_sum[l][p] << " "
               << l << " " << low[0] << " " << low[1] << " " << low[2] << " "
               << high[0] << " " << high[1] << " " << high[2] << std::endl;
      }
    }

    //__________________________________
    //  debugging output
    if (d_myworld->myRank() == 0) {
      // calculate total cost for normalization
      double total_measured               = 0;
      double total_predicted              = 0;
      double total_measured_error         = 0;
      double total_measured_percent_error = 0;
      double total_volume                 = 0;

      for (int r = 0; r < (int)regions.size(); r++) {
        total_measured += measured_sum[l][r];
        total_predicted += predicted_sum[l][r];
        total_measured_error += fabs(measured_sum[l][r] - predicted_sum[l][r]);
        total_measured_percent_error +=
          fabs(measured_sum[l][r] - predicted_sum[l][r]) / measured_sum[l][r];
        total_volume += regions[r].getVolume();
      }

      stats << "Profile Error: " << l << " "
            << total_measured_error / regions.size() << " "
            << total_measured_percent_error / regions.size() << " "
            << total_measured << " " << total_predicted << std::endl;
    }
  }

  smpe /= num_patches;
  smape /= num_patches;

  if (d_myworld->myRank() == 0 && stats.active()) {
    std::cout << "SMPE: " << smpe << " sMAPE: " << smape
              << " MAXsPE: " << max_error << std::endl;
  }

  double meanCostm = 0;
  double meanCostp = 0;
  double maxCostm  = proc_costsm[0];
  double maxCostp  = proc_costsp[0];
  int maxLocm      = 0;
  int maxLocp      = 0;

  for (int p = 0; p < d_myworld->nRanks(); p++) {
    meanCostm += proc_costsm[p];
    meanCostp += proc_costsp[p];

    if (maxCostm < proc_costsm[p]) {
      maxCostm = proc_costsm[p];
      maxLocm  = p;
    }

    if (maxCostp < proc_costsp[p]) {
      maxCostp = proc_costsp[p];
      maxLocp  = p;
    }
  }

  if (d_myworld->myRank() == 0) {
    stats << "LoadBalance Measured:  Mean:" << meanCostm / d_myworld->nRanks()
          << " Max:" << maxCostm << " on processor " << maxLocm << std::endl;
    stats << "LoadBalance Predicted:  Mean:" << meanCostp / d_myworld->nRanks()
          << " Max:" << maxCostp << " on processor " << maxLocp << std::endl;
  }
#if 0
  if(maxCostm/maxCostp>1.1)
  {
     std::stringstream str;
    if(d_myworld->myRank()==0)
    {
      stats << d_myworld->myRank() << " Error measured/predicted do not line up, patch cost processor " << maxLocm << " patches:\n";
    }

    //for each level
    for (int l=0; l<currentGrid->numLevels();l++)
    {
      LevelP level=currentGrid->getLevel(l);
      std::vector<Region> regions(level->numPatches());

      for(int p=0; p<level->numPatches();p++)
      {
        const Patch *patch=level->getPatch(p);
        regions[p]=Region(patch->getCellLowIndex(),patch->getCellHighIndex());
      }

      double maxError=0;
      int maxLoc=0;
      for(unsigned int r=0; r<regions.size();r++)
      {
        const Patch *patch=level->getPatch(r);

        int proc=d_lb->getPatchwiseProcessorAssignment(patch);

        if(proc==maxLocm)
        {
          if(d_myworld->myRank()==0)
          {
            stats << "    level: " << l << " region: " << regions[r] << " measured:" << measured_sum[l][r] << " predicted:" << predicted_sum[l][r] << std::endl;
          }
          //find max error region
          double error=measured_sum[l][r]-predicted_sum[l][r];
          if(error>maxError)
          {
            maxError=error;
            maxLoc=r;
          }
        }
      }

      IntVector low=regions[maxLoc].getLow()/d_minPatchSize[l];
      IntVector high=regions[maxLoc].getHigh()/d_minPatchSize[l];

      if(d_myworld->myRank()==0)
        str << "        map entries for region:" << regions[maxLoc] << std::endl;

      //loop through datapoints
      for(CellIterator iter(low,high); !iter.done(); iter++)
      {
        //search for point in the map
        std::map<IntVector,Contribution>::iterator it=costs[l].find(*iter);

        //if in the map
        if(it!=costs[l].end())
        {
         str << "              " << d_myworld->myRank() << " level: " << l << " key: " << it->first << " measured: " << it->second.current << " predicted: " << it->second.weight << std::endl;
        }
      }
    }
    for(int p=0;p<d_myworld->nRanks();p++)
    {
      MPI_Barrier(d_myworld->getComm());
      if(p==d_myworld->myRank())
        stats << str.str();
    }
  }
#endif
}

//______________________________________________________________________
//
void
ProfileDriver::finalizeContributions(const GridP currentGrid)
{
  // if(d_myworld->myRank()==0)
  //   std::cout << "Finalizing Contributions in cost profiler on timestep: " <<
  //   d_timesteps << "\n";

  if (stats.active()) {
    outputError(currentGrid);
  }

  d_timesteps++;
  // for each level
  for (int l = 0; l < (int)costs.size(); l++) {
    // for each datapoint
    for (std::map<IntVector, Contribution>::iterator iter = costs[l].begin();
         iter != costs[l].end();) {
      // copy the iterator (so it can be deleted if we need to)
      std::map<IntVector, Contribution>::iterator it = iter;

      // create a reference to the data
      Contribution& data = it->second;

      // update timestep
      if (data.current > 0) {
        data.timestep = 0;
      } else {
        data.timestep++; // this keeps track of how long it has been since the
                         // data has been updated on this processor
      }

      if (d_timesteps <= 2) {
        // first couple timesteps should be set to last timestep to initialize
        // the system the first couple timesteps are not representative of the
        // actual cost as extra things are occuring.
        // if(data.current==0)
        //   std::cout << d_myworld->myRank() << "WARNING current is equal to 0
        //   at " << it->first << " weight: " << data.weight << " current: " <<
        //   data.current << std::endl;
        data.weight   = data.current;
        data.timestep = 0;
      } else {
        if (d_type == MEMORY) {
          // update exponential average      See Equation 5.4 of Justin's
          // Dissertation
          data.weight = d_alpha * data.current + (1 - d_alpha) * data.weight;
        } else // TYPE IS KALMAN               See Eqs. 5.6 to 5.9
        {
          double m = data.p + d_phi;
          double k = m / (m + d_r);
          data.p   = (1 - k) * m; // computing covariance

          // cout << setprecision(12);
          // cout << "m: " << m << " k:" << k << " p:" << data.p << std::endl;

          data.weight = data.weight + k * (data.current - data.weight);
        }
      }

      // reset current
      data.current = 0;

      // increment the iterator (doing it here because erasing might invalidate
      // the iterator)
      iter++;

      // if the data is empty or old
      if (data.weight == 0 ||
          data.timestep > log(.001 * d_alpha) / log(1 - d_alpha)) {
        // cout << d_myworld->myRank() << " erasing data on level " << l << " at
        // index: " << it->first << " data.timestep: " << data.timestep << "
        // threshold: " << log(.001*d_alpha)/log(1-d_alpha) << std::endl; erase saved
        // iterator in order to save space and time
        costs[l].erase(it);
      }
    }
  }
}
//______________________________________________________________________
//
void
ProfileDriver::getWeights(int l,
                          const std::vector<Region>& regions,
                          std::vector<double>& weights)
{
  if (regions.size() == 0) {
    return;
  }

  weights.resize(regions.size());
  std::vector<double> partial_weights(regions.size(), 0);

  for (int r = 0; r < (int)regions.size(); r++) {
    // coarsen region by minimum patch size
    IntVector low  = regions[r].getLow() / d_minPatchSize[l];
    IntVector high = regions[r].getHigh() / d_minPatchSize[l];

    // loop through datapoints
    for (CellIterator iter(low, high); !iter.done(); iter++) {
      // search for point in the map so we don't inadvertantly add it to the
      // map.  If it doesn't exist on this processor it should exist on another
      // processor
      std::map<IntVector, Contribution>::iterator it = costs[l].find(*iter);

      // if in the map
      if (it != costs[l].end()) {
        // add cost to weight
        partial_weights[r] += it->second.weight;
      }
    }
  }

  // allreduce sum weights
  if (d_myworld->nRanks() > 1) {
    Uintah::MPI::Allreduce(&partial_weights[0],
                           &weights[0],
                           weights.size(),
                           MPI_DOUBLE,
                           MPI_SUM,
                           d_myworld->getComm());
  }
}

//______________________________________________________________________
//
void
ProfileDriver::initializeWeights(const Grid* oldgrid, const Grid* newgrid)
{
  if (d_timesteps == 0) {
    return;
  }

  //__________________________________
  ////for each level
  for (int l = 1; l < newgrid->numLevels(); l++) {
    std::vector<Region> old_level;

    if (oldgrid && oldgrid->numLevels() > l) {
      // create old_level vector
      for (int p = 0; p < oldgrid->getLevel(l)->numPatches(); p++) {
        const Patch* patch = oldgrid->getLevel(l)->getPatch(p);
        old_level.push_back(
          Region(patch->getCellLowIndex(), patch->getCellHighIndex()));
      }
    }

    // get weights on old_level
    std::vector<double> weights;
    getWeights(l, old_level, weights);

#if 1
    double volume = 0;
    double weight = 0;
    // compute average cost per datapoint
    for (int r = 0; r < (int)old_level.size(); r++) {
      volume += old_level[r].getVolume();
      weight += weights[r];
    }
    double initial_cost = weight / volume * d_minPatchSizeVolume[l];
#elif 1
    double weight = DBL_MAX;
    // compute minimum cost per datapoint
    for (int r = 0; r < (int)old_level.size(); r++) {
      if (weights[r] / old_level[r].getVolume() < weight) {
        weight = weights[r] / old_level[r].getVolume();
      }
    }
    double initial_cost = weight * d_minPatchSizeVolume[l] * 2;
#else
    double weight = 0;
    // compute maximum cost per datapoint
    for (int r = 0; r < (int)old_level.size(); r++) {
      if (weights[r] / old_level[r].getVolume() > weight) {
        weight = weights[r] / old_level[r].getVolume();
      }
    }
    double initial_cost = weight * d_minPatchSizeVolume[l];
#endif
    // if there is no cost data
    if (initial_cost == 0) {
      // set each datapoint to the same small value
      initial_cost = 1;
    }

    // compute regions in new level that are not in old

    std::vector<Region> new_regions_partial, new_regions, dnew,
      dold(old_level.begin(), old_level.end());

    // create dnew to contain a subset of the new patches
    for (int p = 0; p < newgrid->getLevel(l)->numPatches(); p++) {
      const Patch* patch = newgrid->getLevel(l)->getPatch(p);
      // distribute regions accross processors in order to parallelize the
      // differencing operation
      if (p % d_myworld->nRanks() == d_myworld->myRank()) {
        dnew.push_back(
          Region(patch->getCellLowIndex(), patch->getCellHighIndex()));
      }
    }

    // compute difference
    new_regions_partial = Region::difference(dnew, dold);

    // gather new regions onto each processor (needed to erase old data)

    int mysize = new_regions_partial.size();
    std::vector<int> recvs(d_myworld->nRanks(), 0),
      displs(d_myworld->nRanks(), 0);

    // gather new regions counts
    if (d_myworld->nRanks() > 1) {
      Uintah::MPI::Allgather(&mysize,
                             1,
                             MPI_INT,
                             &recvs[0],
                             1,
                             MPI_INT,
                             d_myworld->getComm());
    } else {
      recvs[0] = mysize;
    }

    int size = recvs[0];
    recvs[0] *= sizeof(Region);

    for (int p = 1; p < d_myworld->nRanks(); p++) {
      // compute displacements
      displs[p] = displs[p - 1] + recvs[p - 1];

      // compute number of regions
      size += recvs[p];

      // convert to bytes
      recvs[p] *= sizeof(Region);
    }

    new_regions.resize(size);
    // gather the regions
    if (d_myworld->nRanks() > 1) {
      Uintah::MPI::Allgatherv(&new_regions_partial[0],
                              recvs[d_myworld->myRank()],
                              MPI_BYTE,
                              &new_regions[0],
                              &recvs[0],
                              &displs[0],
                              MPI_BYTE,
                              d_myworld->getComm());
    } else {
      new_regions.swap(new_regions_partial);
    }
#if 0
    if(d_myworld->myRank()==0)
    {
      std::cout << " Old Regions: ";
      for(vector<Region>::iterator it=old_level.begin();it!=old_level.end();it++)
      {
        IntVector low=it->getLow()/d_minPatchSize[l];
        IntVector high=it->getHigh()/d_minPatchSize[l];
        for(CellIterator iter(low,high); !iter.done(); iter++)
          std::cout << *iter << " ";
      }
      std::cout << std::endl;
      std::cout << " New Regions: ";
      for(vector<Region>::iterator it=new_regions.begin();it!=new_regions.end();it++)
      {
        IntVector low=it->getLow()/d_minPatchSize[l];
        IntVector high=it->getHigh()/d_minPatchSize[l];
        for(CellIterator iter(low,high); !iter.done(); iter++)
          std::cout << *iter << " ";
      }
      std::cout << std::endl;
    }
#endif

    int p = 0;
    // initialize weights
    for (std::vector<Region>::iterator it = new_regions.begin();
         it != new_regions.end();
         it++) {
      // add regions to my map
      IntVector low  = it->getLow() / d_minPatchSize[l];
      IntVector high = it->getHigh() / d_minPatchSize[l];

      // loop through datapoints
      for (CellIterator iter(low, high); !iter.done(); iter++) {

        std::map<IntVector, Contribution>::iterator it = costs[l].find(*iter);

        // erase any old data in map
        if (it != costs[l].end()) {
          costs[l].erase(it);
        }

        if (p++ % d_myworld->nRanks() ==
            d_myworld->myRank()) // distribute new regions accross processors
        {
          // add cost to current contribution
          costs[l][*iter].weight = initial_cost;
          //          if(l==2)
          //            std::cout << " initializing " << *iter*d_minPatchSize[l]
          //            << " to " << initial_cost << std::endl;
        }
      } // end cell iteration
    }   // end region iteration
  }     // end levels iteration
}

void
ProfileDriver::reset()
{
  for (int i = 0; i < (int)costs.size(); i++) {
    costs[i].clear();
  }
  d_timesteps = 0;
}
