---
title:  "Parallel domain decomposition for particle methods: Part 4"
subheadline: "Biswajit Banerjee"
description: "Applying the Plimpton method for migrating particles"
date:  2017-07-27 10:30:00
categories:
    - MPI
    - C++
excerpt_separator: <!--more-->
toc: true
toc_label: "Contents"
toc_sticky: true
---

#### Introduction ####
The Plimpton scheme of communicating ghost information between patches was described
in [Part 3]({{ site.baseurl }}/mpi/c++/parallel-domain-decomposition-part-3/) of this series.
Let us now see how a similar approach can be used to migrate particles that have
moved across patches.
<!--more-->

In the animation below we just move the particles within each patch randomly.  To make
the identity of the particles clear, we have used solid circles for the patch particles and
three-quarter circle for the particles in the ghost regions. As you can see, some of the
particles have moved outside the patches and need either to be deleted (if they have
left the computational domain - assuming that the domain size remains unchanged) or they
need to be moved to adjacent patches.

<div align="center" style="border:1px solid black">
<div>
  <input name="restartMotion" type="button" value="Repeat particle motion" onclick="particleMotion.restartAnimation()" />
</div>
<div>
  <canvas id="particle-motion" height="500" width="500"></canvas>
</div>
</div>
<p/>

#### Plimpton's scheme for migrating particles ####
If we run Plimpton's scheme in reverse order, we can move the particles to the appropriate
patches with only four communication steps in 2D and six in 3D.  Notice in the animation
below that we start with a search region in the x-direction that contains the top and bottom
patches along with the right (or left) patch.  We relocate particles in this region first and
then need to move particles only in the top and bottom patches.  Note also that the ghost particles
have been moved back to their original locations, indicating that we can ignore these during the
migration process.  Depending on the requirements of the problem, we may either delete particles
that have left the domain, introduce them back in a periodic manner, or extend the domain itself.

<div align="center" style="border:1px solid black">
<div>
  <input name="restartMigrate" type="button" value="Migrate particles across patches" onclick="particleMigrate.restartAnimation()" />
</div>
<div>
  <canvas id="particle-migrate" height="500" width="500"></canvas>
</div>
</div>
<p/>

#### MPI implementation ####
The implementation of the migration process is similar to that for ghost exchange. A typical
`migrateParticles` function can have the following form:

{% highlight cpp %}
using ParticleIDHashMap = std::unordered_set<ParticleID>;
void migrateParticle(...., const Vec& patchWidth, ParticlePArray& patchParticles)
{
  ParticleIDHashMap sentParticles;  // sent particles per process
  ParticlePArray    recvParticles;  // received particles per process
  // First migrate in the x-direction
  d_patchP->sendRecvMigrateXMinus(boostWorld, patchWidth, patchParticles);
  d_patchP->sendRecvMigrateXPlus(boostWorld, patchWidth, patchParticles);
  d_patchP->waitToFinishX();
  d_patchP->combineSentParticlesX(sentParticles);
  d_patchP->combineReceivedParticlesX(recvParticles);
  d_patchP->deleteSentParticles(sentParticles, patchParticles);
  d_patchP->addReceivedParticles(recvParticles, patchParticles);
  // Next migrate in the y-direction
  sentParticles.clear();
  recvParticles.clear();
  d_patchP->sendRecvMigrateYMinus(boostWorld, patchWidth, patchParticles);
  d_patchP->sendRecvMigrateYPlus(boostWorld, patchWidth, patchParticles);
  d_patchP->waitToFinishY();
  d_patchP->combineSentParticlesY(sentParticles);
  d_patchP->combineReceivedParticlesY(recvParticles);
  d_patchP->deleteSentParticles(sentParticles, patchParticles);
  d_patchP->addReceivedParticles(recvParticles, patchParticles);
  // Next migrate in the z-direction
  sentParticles.clear();
  recvParticles.clear();
  d_patchP->sendRecvMigrateZMinus(boostWorld, patchWidth, patchParticles);
  d_patchP->sendRecvMigrateZPlus(boostWorld, patchWidth, patchParticles);
  d_patchP->waitToFinishZ();
  d_patchP->combineSentParticlesZ(sentParticles);
  d_patchP->combineReceivedParticlesZ(recvParticles);
  d_patchP->deleteSentParticles(sentParticles, patchParticles);
  d_patchP->addReceivedParticles(recvParticles, patchParticles);
  // delete outgoing particles (if needed by the problem)
  d_patchP->removeParticlesOutsidePatch(patchParticles);
}
{% endhighlight %}

In [Part 2]({{ site.baseurl }}/mpi/c++/parallel-domain-decomposition-part-2/) we defined
a `PatchNeighborComm` struct and a `Patch` struct.  We can keep the `PatchNeighborComm`
struct in the same form, with the possible addition of a method of two.  However, the `Patch`
struct becomes considerably simplified as show below.

##### Patch struct #####
The `Patch` struct described in [Part 3]({{ site.baseurl }}/mpi/c++/parallel-domain-decomposition-part-3/)
now has a few more methods.  Let us see how some of these new functions may be implemented.

The first new function is `sendRecvMigrateXMinus` which is the equivalent of `sendRecvGhostXMinus`
for th emigration process.  Note that the only difference between these two function is the
definition of the search box.

{% highlight cpp %}
void 
Patch::sendRecvMigrateXMinus(boost::mpi::communicator& boostWorld, 
                             const Vec& neighborWidth,
                             const ParticlePArray& particles) 
{
  if (d_xMinus.d_boundary == PatchBoundary::inside) {
    Vec neighborLower = d_lower - neighborWidth;
    Vec neighborUpper = d_upper + neighborWidth;
    neighborUpper.setX(d_lower.x());
    Box neighborBox(neighborLower, neighborUpper);
    d_xMinus.asyncSendRecv(boostWorld, 
                           neighborBox, d_tolerance,
                           particles);
  }
}
{% endhighlight %}

The next new method is `combineSentParticlesX` which is defined as:

{% highlight cpp %}
void 
Patch::combineSentParticlesX(ParticleIDHashMap& sent) 
{
  d_xMinus.combineSentParticles(sent);
  d_xPlus.combineSentParticles(sent);
}
{% endhighlight %}

The `combineSentParticles` method in `PatchNeighborComm` is defined as

{% highlight cpp %}
void 
PatchNeighborComm::combineSentParticles(ParticleIDHashMap& sent) 
{
  if (!d_sentParticles.empty()) {
    for (const auto& particle : d_sentParticles) {
      sent.insert(particle->getId());
    }
  }
}
{% endhighlight %}

One also needs to delete the sent particles from the patch, using the method
`deleteSentParticles`;  this is where the use of a hashmap becomes handy.

{% highlight cpp %}
void
Patch::deleteSentParticles(const ParticleIDHashMap& sent,
                           ParticlePArray& particles)
{
  if (sent.size() > 0) {
    particles.erase(
      std::remove_if(
        particles.begin(), particles.end(),
        [&sent](const ParticleP& particle) {
          return (sent.find(particle->getId()) != std::end(sent));
        }),
      std::end(particles));
  }
}
{% endhighlight %}

Finally, we add the received particles to the list of particles in the patch using
`addReceivedParticles`:

{% highlight cpp %}
void
Patch::addReceivedParticles(const ParticlePArray& received,
                            ParticlePArray& particles)
{
  particles.insert(particles.end(), received.begin(), received.end());
}
{% endhighlight %}

In some special cases, we will also need to remove particles outside the domain.  One
approach is to use `removeParticlesOutsidePatch`, but this step is typically not recommended
in general as it is costly and often not necessary.

{% highlight cpp %}
void 
Patch::removeParticlesOutsidePatch(ParticlePArray& particles)
{
  Box box(d_lower, d_upper);
  double epsilon = d_tolerance;
  particles.erase(
    std::remove_if(
      particles.begin(), particles.end(),
      [&box, &epsilon](ParticleP particle) {
        if (box.inside(particle->currentPos(), epsilon)) {
          return false;
        }
        return true;
      }),
    particles.end());
}
{% endhighlight %}
<p/>


#### Remarks ####
In the next part of this series we will explore how information about forces can be communicated
across patches.

<script src="{{ site.baseurl }}/assets/js/d3.v4.min.js"></script>
<script src="{{ site.baseurl }}/assets/js/colorbrewer.min.js"></script>
<script src="{{ site.baseurl }}/assets/js/seedrandom.min.js"></script>
<script src="{{ site.baseurl }}/assets/js/particleMigrate.js"></script>
<script src="{{ site.baseurl }}/assets/js/particleMotion.js"></script>

