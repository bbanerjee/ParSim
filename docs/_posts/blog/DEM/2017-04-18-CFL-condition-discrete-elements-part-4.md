---
title:  "The CFL condition for explicit discrete element methods:4"
subheadline: "Biswajit Banerjee"
description: "Part 4: Stability condition and impact"
date:  2017-04-18 10:30:00
categories:
    - DEM
excerpt_separator: <!--more-->
toc: true
toc_label: "Contents"
toc_sticky: true
---


#### Introduction ####
In [Part 3]({{ site.baseurl }}/dem/CFL-condition-discrete-elements-part-3/) of this article I
discussed the approach where the equations for a system of rigid bodies are
approximated by a spring-mass system.  The numerical stability conditions of that
system are then taken to be representative of the system of discrete rigid bodies.
<!--more-->
In that case, for a typical central difference scheme, the time step size is
<div>
$$
  \Delta t \le \frac{2}{\omega}(\sqrt{\zeta^2+1} - \zeta) 
$$
</div>
In terms of the spring stiffness $$k$$, the damping coefficient $$c$$, and the mass $$m$$,
we have
<div>
$$
  \Delta t \le 2\sqrt{\frac{m}{k}}(\sqrt{\frac{c^2}{4km}+1} - \frac{c}{2\sqrt{km}}) 
$$
</div>
{:.notice}
This condition assumes that the system of ODEs is linear.  Note that this assumption
is clearly violated for complex discrete element calculations.

Let us now look at some simple stability bounds for a system of discontinuous
rigid bodies that interact occasionally via impact.

####  One-dimensional equations for impact ####
The impact between two rigid bodies in discrete elements can be simplified to the following
one-dimensional system of equations:
<div>
$$
  m\ddot{d} = \begin{cases}
                -c\dot{d} - kd & \quad \text{for} \quad d \ge 0 \quad \text{:contact} \\
                0 & \quad \text{for} \quad d \lt 0 \quad \text{:separated} 
              \end{cases}
$$
</div>
where $$d$$ is the gap distance.  An animation of the one-dimensional problem is shown below.
In the animation, the smaller ball has been replaced with a spring and damper system (shown
as a spring in the figure).
Note that while the bodies are in contact, the governing equation is identical to
that for a system of springs-dampers-masses.

<div>
  <canvas id="ballball" width="500" height="300"></canvas>
</div>

At the end of the impact event, the bodies will typically separate.  The duration
of impact (which can be found by solving the linear ODE for $$d \ge 0$$), is
<div>
$$
  t_d = \frac{\pi}{\sqrt{1 - \zeta^2}}
$$
</div>
and the exit velocity is
<div>
$$
  v_e = v_0 \exp\left[-\frac{\zeta\pi}{\sqrt{1-\zeta^2}}\right]
$$
</div>

####  Stability of central difference for one-dimensional impact ####
An algorithm that is used to solve the above non-linear problem will be stable is
the exit velocity satisfies the condition $$v_e \le v_0$$ where $$v_0$$ is the
relative velocity of the two bodies.

Feng's work (IJNME 2005, 64:1959) used numerical solutions of the one-dimensional problem to prove that
discrete element calculations can become unconditionally unstable for
$$\Delta t < 2/\omega$$, except for a few values where the exit velocity is equal
to $$v_0$$.  However, a small amount of damping can improve the stability behavior
of the system.

#### Remarks ####
In general, there is no known closed-form expression that can be used to determine
the region of stability of the nonlinear system of ODEs that describes DEM calculations.
As a result, practitioners typically use timesteps that are a small fraction
of that suggested by the analysis of the continuous linear ODE that describes
a spring-mass-damper system.

<script src="{{ site.baseurl }}/assets/js/d3.v4.min.js"></script>
<script src="{{ site.baseurl }}/assets/js/demImpact.js"></script>
