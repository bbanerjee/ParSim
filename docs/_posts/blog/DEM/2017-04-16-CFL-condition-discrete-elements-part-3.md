---
layout: posts
title:  "The CFL condition for explicit discrete element methods:3"
subheadline: "Biswajit Banerjee"
description: "Part 3: Stability condition for discrete elements"
date:  2017-04-16 10:30:00
categories:
    - DEM
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

#### Introduction ####
In [Part 1]({{ site.url }}/dem/CFL-condition-discrete-elements/) of this article,
we revisited the CFL condition and in [Part 2]({{ site.url }}/dem/CFL-condition-discrete-elements-part-2/) we showed how the CFL condition and the von Neumann stability condition are
identical for linear wave equation PDEs.

Both these conditions apply to PDEs for continuous systems (with the possibility
of shocks) that have been discretized.  However, the discrete element method is
discretizes a set of coupled ODEs that apply to each rigid particle in a system.
Does it even make sense to talk about the CFL condition for discrete elements?

I suggest we should just use the term "stability condition" for discrete element 
methods. Even though a system of discrete particles can exhibit wave-like behavior,
it is not straightforward to connect the behavior of the system to an equivalent
PDE that can be analyzed using the CFL/von Neumann techniques.

First, let us explore the ODEs governing a single rigid body and then examine the
numerical stability of the ODE during integration.

#### Euler's laws of motion ####
Discrete element methods model systems of objects as rigid bodies that interact
via contact.  The motion of a single rigid body with respect to an **inertial frame**
can be described using **Euler's equations of motion**:
<div>
$$
  \mathbf{f} = \frac{d\mathbf{p}}{dt} = m \frac{d\mathbf{v}}{dt} \quad \text{and} \quad
  \boldsymbol{\tau} = \frac{d\mathbf{l}}{dt} = \frac{d}{dt}(\boldsymbol{I}\cdot\boldsymbol{\Omega})
$$
</div>
where $$\mathbf{p}$$ is the linear momentum, $$\mathbf{l}$$ is the angular momentum,
$$m$$ is the mass, $$\boldsymbol{I}$$ is the mass moment of inertia, $$\mathbf{v}$$ is
the linear velocity, $$\boldsymbol{\Omega}$$ is the angular velocity, $$\mathbf{f}$$ is
the force, $$\boldsymbol{\tau}$$ is the torque, and $$t$$ is the time.  Note that
in this case all quantities are computed with respect to the center of mass of the body.

For numerical purposes, it is common to express the linear and angular velocities in
terms of derivatives of spatial position ($$\mathbf{x}$$) and
orientation ($$\mathbf{\theta}$$) as
<div>
$$
  \mathbf{v} = \frac{d\mathbf{x}}{dt} \quad \text{and} \quad
  \boldsymbol{\Omega} = \frac{d\boldsymbol{\theta}}{dt}
$$
</div>

#### Discretization of Euler equations ####
Note that both Euler equations have the form
<div>
$$
  \ddot{\mathbf{x}}(t) = \mathbf{A}(\mathbf{x}(t))
$$
</div>
where $$\ddot{a} := d^2 a/dt^2$$.  If we know the initial conditions $$\mathbf{x}(0) = \mathbf{x}_0$$ and $$\dot{\mathbf{x}}(0) = \mathbf{v}_0$$, we can integrate the two equations
simultaneously using a [Velocity-Verlet](https://en.wikipedia.org/wiki/Verlet_integration)
scheme:
<div>
$$
  \begin{align}
  \mathbf{v}_{n+1/2} & = \mathbf{v}_n + \tfrac{1}{2} \dot{\mathbf{v}}_n \Delta t 
  & \boldsymbol{\Omega}_{n+1/2} & = \boldsymbol{\Omega}_n +
     \tfrac{1}{2} \dot{\boldsymbol{\Omega}}_n \Delta t \\
  \mathbf{x}_{n+1} & = \mathbf{x}_n + \mathbf{v}_{n+1/2} \Delta t
  & \boldsymbol{\theta}_{n+1} & = \boldsymbol{\theta}_n +
    \boldsymbol{\Omega}_{n+1/2} \Delta t \\
  \mathbf{v}_{n+1} & = \mathbf{v}_{n+1/2} + \tfrac{1}{2} \dot{\mathbf{v}}_{n+1} \Delta t 
  & \boldsymbol{\Omega}_{n+1} & = \boldsymbol{\Omega}_n +
     \tfrac{1}{2} \dot{\boldsymbol{\Omega}}_{n+1} \Delta t 
  \end{align}
$$
</div>
Note that the angular momentum equation is typically
expressed in body coordinates to avoid a time-varying moment of inertia and solved
accordingly.
{:.notice--warning}

A more common approach is to just use a standard **central difference** approach, i.e.,
<div>
$$
  \begin{align}
  \mathbf{v}_{n+1/2} & = \mathbf{v}_{n-1/2} + \dot{\mathbf{v}}_n \Delta t 
  & \boldsymbol{\Omega}_{n+1/2} & = \boldsymbol{\Omega}_{n-1/2} +
     \dot{\boldsymbol{\Omega}}_n \Delta t \\
  \mathbf{x}_{n+1} & = \mathbf{x}_n + \mathbf{v}_{n+1/2} \Delta t
  & \boldsymbol{\theta}_{n+1} & = \boldsymbol{\theta}_n +
    \boldsymbol{\Omega}_{n+1/2} \Delta t 
  \end{align}
$$
</div>

#### Stability of central difference scheme ####
It is not very easy to explore the stability of the general second-order ODE in the presence
of arbitrary external forces and torques.  Instead, we can draw upon the large literature
on the subject from contact problems in finite element analysis (see [Nonlinear Finite Elements for Continua and Structures](https://books.google.co.nz/books?id=BQpfAQAAQBAJ&dq=nonlinear+finite+elements+of+continua+and+structures), Chapter 6, for an excellent exposition) and explore the more relevant problem of contact between two rigid bodies.

For simplicity let us follow the approach taken by Y. T. Feng
(*On the central difference algorithm in discrete element modelling of impact*, IJNME 64:1959-1980, 2005).
We assume, as is popular in discrete element methods, that the contact between two rigid bodies
can be modeled using a one-dimensional spring-damper system where the spring has a stiffness $$k$$
and the damper has a damping coefficient $$c$$.

##### One-dimensional unforced system #####
If there are no external forcings, the two-body system can be described by the equation
<div>
$$
  m\ddot{x} + c\dot{x} + kx = 0 \,\quad x(0) = 0, \dot{x}(0) = v_0
$$
</div>
This ODE has solutions of the form
<div>
$$
  x(t) = e^{\alpha t}
$$
</div>
Therefore we have
<div>
$$
  \alpha^2 m + \alpha c + k = 0
  \quad \implies \quad
  \alpha = \frac{-c \pm \sqrt{c^2 - 4 k m}}{2m} 
$$
</div>
If $$c^2 > 4 km$$ the solutions are real and not oscillatory (overdamped).  On the other
hand, if $$c^2 \lt 4 km$$ the solutions are imaginary and therefore oscillatory (underdamped).
The **critical damping** value is given by
<div>
$$
  c_c^2 = 4 k m \quad \equiv \left(\frac{c_c}{2m}\right)^2 = \frac{k}{m} = \omega^2
  \quad \implies \quad c_c = 2\omega m
$$
</div>
where $$\omega$$ is the natural frequency of the system. Then the original ODE can be
written in terms of $$c_c$$ and $$\omega$$ as
<div>
$$
  \ddot{x} + 2\zeta\omega\dot{x} + \omega^2 x = 0 \,\quad \text{where} \quad \zeta := \frac{c}{c_c}
$$
</div>

##### Stability using the Hurwitz matrix approach #####
Using a central difference scheme for second derivatives and an upwind scheme for first derivatives,
we have
<div>
$$
  \frac{x_{n+1} - 2 x_n + x_{n-1}}{\Delta t^2} +
  2\zeta\omega\frac{x_n - x_{n-1}}{\Delta t} + \omega^2 x_n = 0 
$$
</div>
This is a linear difference equation that has solutions of the form $$x_n = r^n$$ (see [Recurrence_relation](https://en.wikipedia.org/wiki/Recurrence_relation) for a proof). Plugging in these solutions, we get
<div>
$$
  \frac{r^{n+1} - 2 r^n + r^{n-1}}{\Delta t^2} +
  2\zeta\omega\frac{r^n - r^{n-1}}{\Delta t} + \omega^2 r^n = 0 
$$
</div>
or
<div>
$$
  \frac{r^2 - 2 r + 1}{\Delta t^2} +
  2\zeta\omega\frac{r - 1}{\Delta t} + \omega^2 r = 0 
$$
</div>
or,
<div>
$$
  r^2 + (g+h-2)r + (1-g) = 0 \quad \text{where} \quad g:= 2\zeta\omega\Delta t ~,~~ h := \omega^2\Delta t^2
$$
</div>
We can now apply the [z-transform](https://en.wikipedia.org/wiki/Bilinear_transform) to this
equation because $$r = \exp(\gamma \Delta t)$$ (see Wikipedia article on recurrence relations)
to get
<div>
$$
 \left(\frac{1+z}{1-z}\right)^2 + (g+h-2)\frac{1+z}{1-z} + (1-g) = 0
$$
</div>
or
<div>
$$
 (1+z)^2 + (g+h-2)(1-z^2) + (1-g)(1-z)^2 = 0
$$
</div>
We can easily generate the [Hurwitz matrix](https://en.wikipedia.org/wiki/Hurwitz_matrix)
($$\mathbf{H}$$) of this
polynomial. Stability requires the leading principal minors of this matrix be positive.
In this case, we have
<div>
$$
  \mathbf{H} = \begin{bmatrix} 2g & 0 \\ 4 - 2g -h & h \end{bmatrix}
    = \begin{bmatrix} 2(2\zeta\omega\Delta t) & 0 \\ 4 - 2(2\zeta\omega\Delta t) -\omega^2\Delta t^2 & \omega^2\Delta t^2 \end{bmatrix}
$$
</div>
and the stability requirements
<div>
$$
  \begin{align}
  g \ge 0  & \implies \quad \zeta\omega\Delta t \ge 0 \\
  4 - 2g -h \ge 0 & \implies \quad 4 - 4\zeta\omega\Delta t - \omega^2\Delta t^2 \ge 0 \\
  h \ge 0 & \implies \quad \omega^2 \Delta t^2 \ge 0
  \end{align}
$$
</div>
For systems with positive damping, the only non-trivial condition above is the second one
which leads to a critical timestep of
<div>
$$
  \Delta t \le \frac{2}{\omega}(\sqrt{\zeta^2+1} - \zeta) \,.
$$
</div>
{:.notice--info}
Variations of this expression are often used to determine the stable time step in DEM
calculations.

#### Remarks ####
In the next part of this series we will explore a slightly more realistic model of
DEM calculations and see what can be said about the stability of that model.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

