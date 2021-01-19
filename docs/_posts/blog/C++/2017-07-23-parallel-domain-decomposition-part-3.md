---
layout: posts
title:  "Parallel domain decomposition for particle methods: Part 3"
subheadline: "Biswajit Banerjee"
description: "The Plimpton method of communicating ghost regions"
date:  2017-07-23 10:30:00
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
In [Part 2]({{ site.url }}/mpi/c++/parallel-domain-decomposition-part-2/) of this series
we showed the direct way of communicating ghost particles between patches.  That approach
requires 26 communication steps per patch in three-dimensions.  In this article we discuss
the approach suggested by Steve Plimpton ("Fast parallel algorithms for short-range molecular 
dynamics", Sandia Report SAND91-1144.UC-405, 1993).

Plimpton's paper has been cited almost 15,000 times since its publication.  Among other
things, the paper explains how the number of communication steps can be reduced to six in
three dimensions.

#### Plimpton's scheme for exchanging particles ####
If you run the animation below, you will notice that the left-right ghost particles are
exchanged first.  The ghost regions in the up-down directions are then extended to include
the left-right ghost regions.  The particles in these enlarged regions are then transferred
in the up-down directions.  Therefore, there are only four communication steps in two-dimensions.
The same process can be used in three-dimensions, leading to only six communication steps.

<div align="center" style="border:1px solid black">
<div>
  <input name="restartExchange" type="button" value="Exchange ghost particles" onclick="particlePlimpton.restartAnimation()" />
</div>
<div>
  <canvas id="particle-exchange-plimpton" height="500" width="500"></canvas>
</div>
</div>
<p/>

#### MPI implementation ####
In [Part 2]({{ site.url }}/mpi/c++/parallel-domain-decomposition-part-2/) we defined
a `PatchNeighborComm` struct and a `Patch` struct.  We can keep the `PatchNeighborComm`
struct in the same form, with the possible addition of a method of two.  However, the `Patch`
struct becomes considerably simplified as show below.

##### Patch struct #####
The `Patch` struct now needs only six `PatchNeighborComm` objects but three `waitToFinish`
and `combineReceivedParticles` methods.

{% highlight cpp %}
struct Patch {
  int d_rank;
  double d_ghostWidth, d_tolerance;
  IntVec d_patchMPICoords;
  Vec d_lower, d_upper;
  PatchNeighborComm d_xMinus, d_yMinus, d_zMinus, d_xPlus, d_yPlus, d_zPlus;
  Patch(MPI_Comm& cartComm,
        int rank, const IntVec& mpiCoords, const Vec& lower, const Vec& upper,
        double ghostWidth, double tolerance);
  void setXMinus(MPI_Comm& cartComm); // A patch boundary is at the domain boundary
  void setXPlus(MPI_Comm& cartComm);
  void setYMinus(MPI_Comm& cartComm);
  // .....
  void setZPlus(MPI_Comm& cartComm);
  // Step 1: Send and receive data from x+ and x-
  void sendRecvGhostXMinus(boost::mpi::communicator& boostWorld,
                           const ParticlePArray& patchParticles);
  void sendRecvGhostXPlus(boost::mpi::communicator& boostWorld,
                          const ParticlePArray& patchParticles);
  void waitToFinishX();
  void combineReceivedParticlesX(ParticlePArray& patchParticles);
  // Step 2: Send and receive data from y+ and y-
  void sendRecvGhostYMinus(boost::mpi::communicator& boostWorld,
                           const ParticlePArray& patchParticles);
  void sendRecvGhostYPlus(boost::mpi::communicator& boostWorld,
                          const ParticlePArray& patchParticles);
  void waitToFinishY();
  void combineReceivedParticlesY(ParticlePArray& patchParticles);
  // Step 3: Send and receive data from z+ and z-
  void sendRecvGhostZMinus(boost::mpi::communicator& boostWorld,
                           const ParticlePArray& patchParticles);
  void sendRecvGhostZPlus(boost::mpi::communicator& boostWorld,
                          const ParticlePArray& patchParticles);
  void waitToFinishZ();
  void combineReceivedParticlesZ(ParticlePArray& patchParticles);
};
{% endhighlight %}

The new implementation of the `Patch` struct is shown below.  The initialization of
the struct is the same as before; as are the `setXMinus` etc. methods.  Also, the
`sendRecvGhostXMinus` and `sendRecvGhostXPlus` methods are the same as before.

{% highlight cpp %}
Patch::Patch(//...) {
  //.....  Same as for direct communication
}
void Patch::setXMinus(MPI_Comm& cartComm) {
  //.....  Same as for direct communication
}
//......
{% endhighlight %}

Next we do the communication of the x+ and x- neighbors and combine the
received particles with the particles in the current patch.

{% highlight cpp %}
// Step 1:  Do the x+ and x- communication
void Patch::sendRecvGhostXMinus(//....) {
  //.....  Same as for direct communication
}
void Patch::sendRecvGhostXPlus(//....) {
  //.....  Same as for direct communication
}
void Patch::waitToFinishX() {
  if (d_xMinus.d_boundary == PatchBoundary::inside) {
    d_xMinus.waitToFinish();
  }
  if (d_xPlus.d_boundary == PatchBoundary::inside) {
    d_xPlus.waitToFinish();
  }
}
void Patch::combineReceivedParticlesX(ParticlePArray& patchParticles) {
  received.clear();
  d_xMinus.combineReceivedParticles(patchParticles);
  d_xPlus.combineReceivedParticles(patchParticles);
}
{% endhighlight %}

Now that the `patchParticles` vector has been updated, we can repeat the
process for the y+ and y- directions.  Note that the size of the ghost
region has been expanded in the negative and positive x-direction.

{% highlight cpp %}
// Step 2:  Do the y+ and y- communication
void Patch::sendRecvGhostYMinus(boost::mpi::communicator& boostWorld, const ParticlePArray& patchParticles)
{
  if (d_yMinus.d_boundary == PatchBoundary::inside) { 
    Vec ghostLower = d_lower - Vec(d_ghostWidth, 0.0, 0.0);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, 0.0, 0.0);
    ghostUpper.setY(d_lower.y() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_yMinus.asyncSendRecv(boostWorld, ghostBox, d_tolerance, patchParticles);
  }
}
void Patch::sendRecvGhostYPlus(boost::mpi::communicator& boostWorld, const ParticlePArray& patchParticles)
{
  if (d_yPlus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, 0.0, 0.0);
    ghostLower.setY(d_upper.y() - d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, 0.0, 0.0);
    Box ghostBox(ghostLower, ghostUpper);
    d_yPlus.asyncSendRecv(boostWorld, ghostBox, d_tolerance, patchParticles);
  }
}
void Patch::waitToFinishY() {
  if (d_yMinus.d_boundary == PatchBoundary::inside) {
    d_yMinus.waitToFinish();
  }
  if (d_yPlus.d_boundary == PatchBoundary::inside) {
    d_yPlus.waitToFinish();
  }
}
void Patch::combineReceivedParticlesY(ParticlePArray& patchParticles) {
  received.clear();
  d_yMinus.combineReceivedParticles(patchParticles);
  d_yPlus.combineReceivedParticles(patchParticles);
}
{% endhighlight %}

Finally, we do the third stage of communication in the z-direction.  The ghost-regions
have now been expanded to contain for the x-, x+ and y-, y+ extensions.
{% highlight cpp %}
// Step 2:  Do the z+ and z- communication
void Patch::sendRecvGhostZMinus(boost::mpi::communicator& boostWorld, const ParticlePArray& patchParticles)
{
  if (d_zMinus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, d_ghostWidth, 0.0);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, d_ghostWidth, 0.0);
    ghostUpper.setZ(d_lower.z() + d_ghostWidth);
    Box ghostBox(ghostLower, ghostUpper);
    d_zMinus.asyncSendRecv(boostWorld, ghostBox, d_tolerance, patchParticles);
  }
}
void Patch::sendRecvGhostZPlus(boost::mpi::communicator& boostWorld, const ParticlePArray& patchParticles)
{
  if (d_zPlus.d_boundary == PatchBoundary::inside) {
    Vec ghostLower = d_lower - Vec(d_ghostWidth, d_ghostWidth, 0.0);
    ghostLower.setZ(d_upper.z() - d_ghostWidth);
    Vec ghostUpper = d_upper + Vec(d_ghostWidth, d_ghostWidth, 0.0);
    Box ghostBox(ghostLower, ghostUpper);
    d_zPlus.asyncSendRecv(boostWorld, ghostBox, d_tolerance, patchParticles);
  }
}
void Patch::waitToFinish() {
  if (d_zMinus.d_boundary == PatchBoundary::inside) {
    d_zMinus.waitToFinish();
  }
  if (d_zPlus.d_boundary == PatchBoundary::inside) {
    d_zPlus.waitToFinish();
  }
}
void Patch::combineReceivedParticlesY(ParticlePArray& patchParticles) {
  received.clear();
  d_zMinus.combineReceivedParticles(patchParticles);
  d_zPlus.combineReceivedParticles(patchParticles);
}
{% endhighlight %}
<p/>

##### The particle exchange function #####
The Plimpton particle exchange function the main simulation can then be simplified to
the following.

{% highlight cpp %}
void ParticleCode::exchangeGhostParticles() {
  // Initialize list of all particles in patch + ghost
  mergeParticleVec.clear();
  mergeParticleVec = particleVec;
  // First the x+/x- directions
  d_patchP->sendRecvGhostXMinus(boostWorld, mergeParticleVec);
  d_patchP->sendRecvGhostXPlus(boostWorld, mergeParticleVec);
  d_patchP->waitToFinishX();
  d_patchP->combineReceivedParticlesX(mergeParticleVec);
  // Next the y+/y- directions
  d_patchP->sendRecvGhostYMinus(boostWorld, mergeParticleVec);
  d_patchP->sendRecvGhostYPlus(boostWorld, mergeParticleVec);
  d_patchP->waitToFinishX();
  d_patchP->combineReceivedParticlesX(mergeParticleVec);
  // Next the z+/z- directions
  d_patchP->sendRecvGhostZMinus(boostWorld, mergeParticleVec);
  d_patchP->sendRecvGhostZPlus(boostWorld, mergeParticleVec);
  d_patchP->waitToFinishX();
  d_patchP->combineReceivedParticlesX(mergeParticleVec);
}
{% endhighlight %}

In this case the number of communication steps is much smaller.  However, there is a wait period
at the end of each communication step that may reduce the benefits of the Plimpton approach
in some situations where the particles are unevenly distributed.

#### Remarks ####
Plimpton's scheme is attractive for its simplicity in communicating ghost particle positions.
However, there are two more important communication steps that need to be considered - the
computation of interparticle forces and the migration of particles between patches.  In the next
part of this series, we will discuss the migration of particles when we use the Plimpton scheme.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

<script src="{{ site.url }}/assets/js/d3.v4.min.js"></script>
<script src="{{ site.url }}/assets/js/colorbrewer.min.js"></script>
<script src="{{ site.url }}/assets/js/particlePlimpton.js"></script>

