---
layout: posts
title:  "Parallel domain decomposition for particle methods: Part 1"
subheadline: "Biswajit Banerjee"
description: "The scatter operation"
date:  2017-07-20 10:30:00
categories:
    - MPI
    - C++
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

#### Introduction ####
For parallel particle codes that have to be written quickly (while retaining flexibility), the
task-based parallelism approach doesn't always work well.  The usual approach that is taken
in those situations is some sort of domain decomposition and a lot of associated fine-grained
code for communication between processes.  One tries to strike the appropriate balance between
communication and computation while making sure that the computation is load-balanced.  As a
rule of thumb, less communication is better.

One approach (among many) in particle-based codes that are being parallelized starting from
a serial version is:

* The creation of the particles on the root/master processor.
* Scattering the particles to various processes.
* Communicating ghost regions at processor boundaries.
* Migrating particles that have crossed processor boundaries to the appropriate process.

In the interest of simplicity, we ignore the communication of interparticle forces.

#### Creating and scattering particles ####
Particles are created on the master process (P0) and then transferred to other parallel
processes during the "scatter" operation.  In the animation below, we assume that there
are nine processes - P0 through P8.  The domain is decomposed into nine squares and the
contents of each square are sent to the appropriate process.

<div>
  <input name="restartButton" type="button" value="Scatter particles" onclick="restartAnimation()" />
</div>
<div>
  <canvas id="particle-scatter" height="500" width="500"></canvas>
</div>

#### MPI implementation ####
A possible MPI implementation of the scattering process is described below.  For convenience
we use the `boost::mpi` wrappers around MPI calls is most cases.  However, some MPI calls
do not have associated Boost calls and we have to use the MPI calls directly.

##### MPI setup #####
The first step is to set up the MPI communicator and determine the rank (and MPI coordinates
in a virtual Cartesian topology) of the current process:

{% highlight cpp %}
#include <boost/mpi.hpp>
........
........

void mpiSetup()
{
  // The Boost MPI communicator
  boost::mpi::communicator boostWorld;
  // The standard MPI communicator
  MPI_Comm mpiWorld = MPI_Comm(boostWorld);
  // Create a 3D Cartesian virtual process topology (3 x 3 x 1)
  int dimensions = 3;
  IntVec mpiProcs = {3, 3, 1};
  IntVec periods = {0, 0, 0};
  int reordering = 0;
  MPI_Comm cartesianComm;
  MPI_Cart_create(mpiWorld, dimensions, mpiProcs.data(), periods.data(), reordering, &cartesianComm);
  // Find the process rank
  int mpiRank;
  MPI_Comm_rank(cartesianComm, &mpiRank);
  // Find the number of processes associated with the communicator
  int mpiSize;
  MPI_Comm_size(cartesianComm, &mpiSize);
  // Find the MPI coordinates of the process
  IntVec mpiCoords;
  MPI_Cart_coords(cartesianComm, mpiRank, dimensions, mpiCoords.data());
  // Save the communicator, rank, coordinates, size etc.
  ........
  ........
}
{% endhighlight %}
In the above, the `IntVec` class is an `std::array<int, 3>`.

##### The scatter operation #####
In the scatter operation, the particles are assigned to each patch and then
sent to the appropriate patches using the asynchronous `isend` operation:
{% highlight cpp %}
/**
 * boostComm     : The boost MPI communicator
 * cartesianComm : The MPI Cartesian communicator
 * dimensions    : The number of dimensions in the virtual topology
 * mpiRank       : Rank of the current process
 * mpiSize       : Number of MPI processes
 * mpiProcs      : Vector containing the number of processes in each dimension
 * domainMin/Max : Minimum/maximum corners of box representing physical domain
 * particles     : Vector of shared pointers to particles
 */
void scatterParticles(boost::mpi::communicator& boostComm,
                      MPI_Comm& cartesianComm,
                      int dimensions, int mpiRank, int mpiSize,
                      const IntVec& mpiProcs,
                      const Vec& domainMin, const Vec& domainMax,
                      ParticlePArray& particles)
{
  // Find the physical dimensions of each patch
  Vec patchWidth = (domainMax - domainMin)/mpiProcs;

  // For the root process
  if (mpiRank == 0) {
    // Create a set of send requests
    boost::mpi::request requests[mpiSize-1];
    // Loop through the number of processors in reverse order
    // so that the root processor rank is accessed last
    ParticlePArray insideParticles;
    for (int rank = mpiSize - 1; rank >= 0; --rank) {
      insideParticles.clear();
      // Find the MPI coordinates of the processor
      IntVec mpiCoords;
      MPI_Cart_coords(cartesianComm, rank, dimensions, mpiCoords.data());
      // Convert these MPI coordinates into physical patch coordinates
      Vec patchLower = domainMin + patchWidth*mpiCoords; 
      Vec patchUpper = patchLower + patchWidth; 
      // Find which particles are contained in the current patch
      findParticles(patchLower, patchUpper, particles, insideParticles);
      if (rank > 0) {
        // Send the particles inside the patch to the appropriate process
        // (asynchronous)
        requests[rank-1] = boostComm.isend(rank, 0, insideParticles);
      } else {
        // All that remains in the root patch are particles that have
        // not been sent to other processors. We just copy the insideParticles
        // to particles
        .....
      }
    }
    // Now wait until all the asynchronous data transfer is complete
    boost::mpi::wait_all(requests, requests + mpiSize - 1);
  } else {
    // Receive data from the root patch
    boostComm.recv(0, 0, particles);
  }
}
{% endhighlight %}
Here `ParticlePArray` is a `std::vector<ParticleP>` and `ParticleP` is
a `std::shared_ptr<Particle>`.  The `Particle` class contains particle
data.  For simplicity, we do not consider the performance implications
of an array of structures (as used in this implementation) versus
a structure of arrays (which is more efficient).

#### Remarks ####
In the next part of this series, we will discuss two approaches for inter-patch communication
for particle-based simulations.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

<script src="{{ site.url }}/assets/js/d3.v4.min.js"></script>
<script src="{{ site.url }}/assets/js/colorbrewer.min.js"></script>
<script src="{{ site.url }}/assets/js/particleScatter.js"></script>

