---
title:  "Parallel domain decomposition for particle methods: Part 2"
subheadline: "Biswajit Banerjee"
description: "Communicating ghost regions the direct way"
date:  2017-07-22 10:30:00
categories:
    - MPI
    - C++
excerpt_separator: <!--more-->
toc: true
toc_label: "Contents"
toc_sticky: true
---

#### Introduction ####
The [previous article]({{ site.baseurl }}/mpi/c++/parallel-domain-decomposition-part-1/) in this series
discussed the scatter operation for moving particles to various processes.  In this second part
of the series we will discuss a commonly used method of communicating information between
processes.  Each process is logically mapped to a "patch".
<!--more-->

In the animation below, particles are generated in the red patch "P0" and then scattered
to the other eight patches.  During a particle-based simulation, some information has to
be transferred between adjacent processes.  The amount of information that has to be
communicated depends on a characteristic length scale that is determined by the particle
algorithm.  In the animation, this length is shown by the "ghost" regions outlined in
a darker shade with dashed borders.

<div align="center" style="border:1px solid black">
<div>
  <input name="restartScatter" type="button" value="Scatter particles" onclick="particleScatterGhost.restartAnimation()" />
</div>
<div>
  <canvas id="particle-scatter-ghost" height="500" width="500"></canvas>
</div>
</div>
<p/>

#### Exchanging particles between processes ####
After the particles have been scattered and ghost regions identified, the particles in
the ghost regions are exchanged as depicted in the animation below.  Notice that, in
addition to the exchange between the left-right and top-bottom patches, the information
at corners of patches also have to communicated to the three adjacent patches for a total
of 8 communication steps. For three-dimensional simulations, 26 such communication steps
are needed for each patch. Also notice that all we are doing is increasing the size of
each patch and including regions of overlap between patches.

<div align="center" style="border:1px solid black">
<div>
  <input name="restartExchange" type="button" value="Exchange ghost particles" onclick="particleExchange.restartAnimation()" />
</div>
<div>
  <canvas id="particle-exchange-ghost" height="500" width="500"></canvas>
</div>
</div>
<p/>

#### MPI implementation ####
To keep things manageable, we create a `PatchNeighborComm` struct for communications
between neighbor patches. We also define a `Patch` struct that takes care of the details for each patch.

##### PatchNeighborComm struct #####
The neighbor communication methods are defined as:

{% highlight cpp %}
enum class PatchBoundary : char {
  xminus, xplus, yminus, yplus, zminus, zplus, inside
};
struct PatchNeighborComm {
  PatchBoundary d_boundary;   // Whether the patch has a neighbor
  int d_rank;                 // Rank of the neighbor
  int d_mpiTag = 0;
  boost::mpi::request d_sendRecvReq[2];
  ParticlePArray d_sentParticles; // For sends to neighbor
  ParticlePArray d_recvParticles; // For receives from neighbor
  void setNeighbor(MPI_Comm& cartComm, int myRank,
                   const IntVec& neighborCoords,
                   PatchBoundary boundaryFlag);
  void asyncSendRecv(boost::mpi::communicator& boostWorld,
                     int myRank, const Box& box, const double& tolerance,
                     const ParticlePArray& particles);
  void findParticlesInBox(const Box& box,
                          const ParticlePArray& particles,
                          const double& tolerance,
                          ParticlePArray& inside);
  void waitToFinish(int myRank);
  void combineReceivedParticles(int myRank, ParticlePArray& received);
};
{% endhighlight %}

An implementation of the functions in this struct is shown below:

{% highlight cpp %}
void PatchNeighborComm::setNeighbor(//....) {
  int neighborRank = -1;
  MPI_Cart_rank(cartComm, neighborCoords.data(), &neighborRank);
  d_rank = neighborRank;
  if (neighborRank > -1) {
    d_boundary = PatchBoundary::inside;
  } else {
    d_boundary = boundaryFlag;
  }
}
void PatchNeighborComm::asyncSendRecv(//...) {
  // Find the particles in the ghost box
  findParticlesInBox(box, particles, tolerance, d_sentParticles);
  // Asynchronous send
  d_sendRecvReq[0] = boostWorld.isend(d_rank, d_mpiTag, d_sentParticles);
  // Immediate asynchronous receive
  d_sendRecvReq[1] = boostWorld.irecv(d_rank, d_mpiTag, d_recvParticles);
}
void PatchNeighborComm::findParticlesInBox(//...) { // Straightforward }
void PatchNeighborComm::waitToFinish(//...) {
  // Wait from processes to receive ghost data
  boost::mpi::wait_all(d_sendRecvReq, d_sendRecvReq + 2);
}
void PatchNeighborComm::combineReceivedParticles(//...) {
  received.insert(received.end(), d_recvParticles.begin(), d_recvParticles.end());
}
{% endhighlight %}


##### Patch struct #####
The `Patch` struct takes care of all the communication needs of each patch.  The
definition I cobbled together is listed below.

{% highlight cpp %}
struct Patch {
  int d_rank;
  double d_ghostWidth, d_tolerance;
  IntVec d_patchMPICoords;
  Vec d_lower, d_upper;
  PatchNeighborComm d_xMinus, d_yMinus, d_zMinus, d_xPlus, d_yPlus, d_zPlus;
  PatchNeighborComm d_xMinus_yMinus, d_xMinus_yPlus, d_xPlus_yMinus, d_xPlus_yPlus;
  PatchNeighborComm d_xMinus_zMinus, d_xMinus_zPlus, d_xPlus_zMinus, d_xPlus_zPlus;
  PatchNeighborComm d_yMinus_zMinus, d_yMinus_zPlus, d_yPlus_zMinus, d_yPlus_zPlus;
  PatchNeighborComm d_xMinus_yMinus_zMinus, d_xMinus_yPlus_zMinus, d_xPlus_yMinus_zMinus, d_xPlus_yPlus_zMinus;
  PatchNeighborComm d_xMinus_yMinus_zPlus, d_xMinus_yPlus_zPlus, d_xPlus_yMinus_zPlus, d_xPlus_yPlus_zPlus;
  Patch(MPI_Comm& cartComm,
        int rank, const IntVec& mpiCoords, const Vec& lower, const Vec& upper,
        double ghostWidth, double tolerance);
  void setXMinus(MPI_Comm& cartComm); // A patch boundary is at the domain boundary
  void setXPlus(MPI_Comm& cartComm);
  void setYMinus(MPI_Comm& cartComm);
  // .....
  void setZPlus(MPI_Comm& cartComm);
  void sendRecvGhostXMinus(boost::mpi::communicator& boostWorld,
                           const ParticlePArray& particles);
  void sendRecvGhostXPlus(boost::mpi::communicator& boostWorld,
                          const ParticlePArray& particles);
  void sendRecvGhostYMinus(boost::mpi::communicator& boostWorld,
                           const ParticlePArray& particles);
  // .....
  void sendRecvGhostZPlus(boost::mpi::communicator& boostWorld,
                          const ParticlePArray& particles);
  // .....
  void sendRecvGhostXMinusYMinus(boost::mpi::communicator& boostWorld,
                                 const ParticlePArray& particles);
  // .....
  void sendRecvGhostXMinusYMinusZminus(boost::mpi::communicator& boostWorld,
                                       const ParticlePArray& particles);
  // .....
  void sendRecvGhostXPlusYPlusZPlus(boost::mpi::communicator& boostWorld,
                                    const ParticlePArray& particles);
  void waitToFinish();
  void combineReceivedParticles(ParticlePArray& received);
};
{% endhighlight %}

The implementation of the `Patch` struct that I came up with is summarized below.
The design can definitely be improved; but recall that our goal is to do a quick
parallelization of an existing serial code.

{% highlight cpp %}
Patch::Patch(//...) {
 d_rank = rank;
 //.....
 setXPlus(cartComm);
 //.....
 setZPlus(cartComm);
}
void Patch::setXMinus(MPI_Comm& cartComm) {
  IntVec neighborCoords = d_patchMPICoords;
  --neighborCoords.x();
  d_xMinus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::xminus);
}
//......
void Patch::setZPlus(MPI_Comm& cartComm) {
  IntVec neighborCoords = d_patchMPICoords;
  ++neighborCoords.z();
  d_zPlus.setNeighbor(cartComm, d_rank, neighborCoords, PatchBoundary::zplus);
}
void Patch::sendRecvGhostXMinus(//....) {
  if (d_xMinus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower;
    Vec ghostUpper = d_upper;
    ghostUpper.setX(d_lower.x() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_xMinus.asyncSendRecv(boostWorld, d_rank, ghostBox, d_tolerance, particles);
  }
}
//......
//......
void Patch::sendRecvGhostXPlusYPlusZPlus(//....) {
//......
}
void Patch::waitToFinish() {
  if (d_xMinus.d_boundary == PatchBoundary::inside) {
    d_xMinus.waitToFinish(d_rank, iteration);
  }
  if (d_xPlus.d_boundary == PatchBoundary::inside) {
    d_xPlus.waitToFinish(d_rank, iteration);
  }
  //.....
}
void Patch::combineReceivedParticles(ParticlePArray& received) {
  received.clear();
  d_xMinus.combineReceivedParticles(d_rank, iteration, received);
  d_xPlus.combineReceivedParticles(d_rank, iteration, received);
  //.....
}
{% endhighlight %}
<p/>

##### The particle exchange function #####
The particle exchange function the main simulation code can then be written as follows.
Note that this design follows the approach taken by Dr. B. Yan for his parallel DEM code
developed a UC Boulder.

{% highlight cpp %}
void
ParticleCode::exchangeGhostParticles() {
  d_patchP->sendRecvGhostXMinus(boostWorld, particleVec);
  d_patchP->sendRecvGhostXPlus(boostWorld, particleVec);
  d_patchP->sendRecvGhostYMinus(boostWorld, particleVec);
  d_patchP->sendRecvGhostYPlus(boostWorld, particleVec);
  d_patchP->sendRecvGhostZMinus(boostWorld, particleVec);
  d_patchP->sendRecvGhostZPlus(boostWorld, particleVec);
  //.....
  d_patchP->sendRecvGhostXPlusYPlusZPlus(boostWorld, particleVec);
  d_patchP->waitToFinish();
  d_patchP->combineReceivedParticles(recvParticleVec);
  mergeParticleVec.clear();
  mergeParticleVec = particleVec;
  mergeParticleVec.insert(mergeParticleVec.end(), recvParticleVec.begin(),
                          recvParticleVec.end());
}
{% endhighlight %}

Clearly, a lot of communication and book-keeping is needed if we follow this approach.  An
alternative approach that uses fewer communication steps is the procedure developed by
Steve Plimpton ("Fast parallel algorithms for short-range molecular dynamics", Sandia Report
SAND91-1144.UC-405, 1993).

#### Remarks ####
In the next part of this series, we will discuss Plimpton's approach for domain decomposition.

<script src="{{ site.baseurl }}/assets/js/d3.v4.min.js"></script>
<script src="{{ site.baseurl }}/assets/js/colorbrewer.min.js"></script>
<script src="{{ site.baseurl }}/assets/js/particleScatterGhost.js"></script>
<script src="{{ site.baseurl }}/assets/js/particleExchange.js"></script>

