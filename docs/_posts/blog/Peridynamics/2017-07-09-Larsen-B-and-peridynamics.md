---
layout: posts
title:  "Can the Larsen-C ice shelf failure be predicted with Peridynamics?"
subheadline: "Biswajit Banerjee"
description: "How good are we at fracture prediction?"
date:  2017-07-09 10:30:00
categories:
    - Fracture
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

#### Introduction ####
On July 7, 2017, the [Project Midas group](http://www.projectmidas.org/blog/multiple-branches/)
released a couple of plots of material velocities and interferograms showing the evolution of a
large rift the [Larsen-C ice shelf](https://en.wikipedia.org/wiki/Larsen_Ice_Shelf#Larsen_C) in Antarctica.  The crack grew 11 miles in a few days before
slowing down.

<div align="center">
<img style="width:400px" alt="Larsen-C cracks" src="{{site.baseurl}}/assets/blogimg/LarsenC_crack.png"/> 
</div>

The image above (from http://www.projectmidas.org/blog/multiple-branches/) shows a number of branches
in the main crack as it approaches the sea.  This image reminded me of crack branching that
we observed in two-dimensional Peridynamics simulations.

<div align="center">
<img style="width:500px" alt="Peridynamics crack" src="{{site.baseurl}}/assets/blogimg/CrackPeri02.png"/> 
</div>

Will it be possible to predict when the main Larsen-C crack will reach the sea using Peridynamics?

####  Issues to be resolved before a simulation ####
I haven't done the simulation yet, but several simplifications are needed before we can
estimate of the time to failure.  Complications arise because:

1. The primary problem is the size of the body that has to be simulated.  The ice-shelf is around
100 km long, 50 km wide, and 0.5 km thick and is supported by sea water.  It is possible that
inertial effects on crack growth are significant in a body of this size.

2. The crack growth rate has varied significantly over time, but is much slower than typical
cracks that are simulated with peridynamics.

3. The crack has reached a stage where branching is prominent.  That implies that the stress state
is much more complex that was the case when the crack appeared to be a tension dominated one.

4.  The geometry of the crack is not linear and the thickness of the ice shelf is not uniform.

5.  Material properties of the ice are not known very well.  For complexities of ice behavior see
[THE MECHANICAL PROPERTIES OF ICE](http://www.dtic.mil/dtic/tr/fulltext/u2/284777.pdf).

Some of these challenges can be addressed by simplifying the problem.  For instance:

1.  We can assume that the body is a cantilevered plate supported by an elastic foundation.

2.  Scale similarity can be used to reduce the size of the problem and also scale the mass
density of the ice to get a simulation that can be completed in a reasonable amount of time.

3.  We can assume isotropic linear elasticity and a constant fracture toughness.

Even with these simplifications, predicting an estimate of the time the crack will take to reach the
sea is not straightforward and the error bounds will be large.

#### Remarks ####
Even though we have seen enormous improvements in algorithms and computational capabilities for
fracture simulation over the last 15 years, predicting dynamic fracture continues to be elusive.
The research that we do at Parresia attempts to improve the predictive fidelity of
numerical methods to solve these types of problems.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

