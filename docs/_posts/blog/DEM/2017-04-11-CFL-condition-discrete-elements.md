---
layout: posts
title:  "The CFL condition for explicit discrete element methods:1"
subheadline: "Biswajit Banerjee"
description: "Part 1: How to estimate a stable timestep size"
date:  2017-04-11 10:30:00
categories:
    - DEM
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---
<link rel='stylesheet' type='text/css' href='{{ site.url }}/assets/js/animateCFL.css' />
- Contents
{:toc}
{:.notice--content}

#### Introduction ####
Solutions of hyperbolic partial differential equations using explicit numerical methods
need a means of limiting the timestep so that the solution is stable.  A criterion 
that is usually used to constrain the step size is the 
[Courant–Friedrichs–Lewy condition](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition).

Let us first explore what the 1928 paper by the Courant, Friedrichs, and Lewy had to say on
stability.  I will then examine the von Neumann form of the condition that is widely used
in forward difference calculations.  Finally, I will say something about the application of these
ideas to discrete element methods.  Depending on my patience, I'll probably break up the article
into two or three parts.

#### The CFL condition ####
The easiest starting point for the analysis is the one-dimensional wave equation:
<div>
$$
  \frac{\partial^2 u}{\partial t^2} - \frac{\partial^2 u}{\partial x^2} = 0 \,.
$$
</div>
Let us assume that the value of the solution $$u(x,t)$$ and its derivatives are known
at time $$t = 0$$.

If we use a grid size $$\Delta x$$ and a timestep $$\Delta t$$, the derivatives can be expressed as
<div>
$$
 \begin{align}
  \frac{\partial^2 u}{\partial x^2} \approx u_{xx} & =
    \frac{1}{\Delta x^2}[u(x+\Delta x, t) - 2u(x, t) + u(x-\Delta x, t)] \\
  \frac{\partial^2 u}{\partial t^2} \approx u_{tt} & =
    \frac{1}{\Delta t^2}[u(x, t+\Delta t) - 2u(x, t) + u(x, t-\Delta t)] \,.
 \end{align}
$$
</div>

Consider the case where $$\Delta t = h$$ and $$\Delta x = kh$$. In
that case, the 1D wave equation can be discretized as
<div>
$$
  u(x, t+h) + u(x, t-h) - u(x+kh, t) + u(x-kh, t) = 0  \,.
$$
</div>

This implies that each point $$x$$ has an evolving region of influence that grows
as shown (for $$k=1$$) in the animation below, i.e., the value at that point has an immediate effect only on
points in the region of influence (also called **domain of dependence**).

<div class="cfl-wave-animation">
</div>

Similarly, there is a **mathematical domain of dependence** of the underlying hyperbolic
PDE that is exact.  What the CFL paper found was that

> The domain of dependence for the difference equation for this mesh will now either lie entirely within
the domain of dependence of the differential equation, or on the other hand will contain this latter
region inside its own domain according as k \< 1 or k \> 1 respectively

This finding implies that if k \< 1, the domain of dependence (DoD) of the discrete equation is inside
the DoD of the original PDE.  You can see that in the animation below (use the slider to
change the value of $$k$$ to 0.5).  In the animation we have kept $$h = \Delta t$$ constant and the
green lines are a proxy for the mathematical domain of dependence.

<div class="cfl-domain-animation">
</div>

If we decrease $$h$$, and keep k \< 1, the domain of dependence of the discrete equation
becomes smaller relative to the mathematical (exact) domain!  Therefore the discrete solution 
will not converge to the exact solution if k \< 1, however small we choose the value of $$h$$ to be.


The CFL paper also proved that convergence is indeed achieved for $$k > 1$$,
or, in this case,
<div>
$$
  \frac{\Delta x}{\Delta t} > 1\,.
$$
</div>

The same observation is true for many different types or linear and nonlinear PDEs.  With that
is mind, We can write the CFL condition as described in the Courant-Friedrich-Lewy paper as follows.

**The CFL Condition**  For each point in the domain of a PDE, the CFL condition is
satisfied if the mathematical domain of dependence is contained in the
numerical domain of dependence.
{:.notice--info}

**The CFL Theorem** The CFL condition is a necessary condition for the **convergence** of a
numerical approximation of a PDE.
The [Lax Equivalence Theorem](https://en.wikipedia.org/wiki/Lax_equivalence_theorem)
then implies that the
CFL condition is also a necessary condition for the **stability** of the numerical scheme if the PDE is linear.
{:.notice--info}

#### Remarks ####
While this is a powerful general result, the constraint is not tight enough for practical 
calculations.  The von Neumann result that will be discussed in the next part of this series
will show how we can come up with a easier way of computing stability conditions for linear PDEs.

The animations are a bit buggy at this point but will do for the purposes of this article.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

<script src="{{ site.url }}/assets/js/d3.v4.min.js"></script>
<script src="{{ site.url }}/assets/js/animateWave.js"></script>
<script src="{{ site.url }}/assets/js/animateCFL.js"></script>

