---
title:  "The CFL condition for explicit discrete element methods:2"
subheadline: "Biswajit Banerjee"
description: "Part 2: The von Neumann version of the CFL condition"
date:  2017-04-14 10:30:00
categories:
    - DEM
excerpt_separator: <!--more-->
toc: true
toc_label: "Contents"
toc_sticky: true
---

<link rel='stylesheet' type='text/css' href='{{ site.baseurl }}/assets/js/animateCFL.css' />

#### Introduction ####
In the [first part]({{ site.baseurl }}/dem/CFL-condition-discrete-elements/) of this article,
we looked at the one-dimensional linear second-order wave equation
<!--more-->
<div>
$$
  \frac{\partial^2 u}{\partial t^2} - \frac{\partial^2 u}{\partial x^2} = 0 
$$
</div>
We saw that for a discretization with grid size $$\Delta x$$ and  timestep $$\Delta t$$, a
central difference approximation is convergent if the CFL condition
<div>
$$
  \frac{\Delta x}{\Delta t} > 1
$$
</div>
is satisfied.
The numerical method is also stable in this case because the PDE is linear.

Of course, in the above equation we have chosen the units of space and time such that the
wave speed is 1.  In general, it is more convenient to work with quantities for which the
wave speed is not 1, and the 1D wave equation has the form
<div>
$$
  \frac{\partial^2 u}{\partial t^2} - c^2\frac{\partial^2 u}{\partial x^2} = 0 \,\quad
  c = \sqrt{\frac{K}{\rho}}
$$
</div>
where, for linear elasticity, $$K$$ is the bulk modulus and $$\rho$$ is the mass-density.
In that case the CFL condition will be scaled accordingly by the wave speed.

Note that the above equation is equivalent to the following first-order wave equations:
<div>
$$
  \frac{\partial u}{\partial t} - c\frac{\partial u}{\partial x} = 0 \,\quad
  \frac{\partial u}{\partial t} + c\frac{\partial u}{\partial x} = 0 
$$
</div>

Linear PDEs with constant coefficients of this form have **plane wave** solutions:
<div>
$$
  u(x, t) = \exp[i(\kappa x - \omega t)]
$$
</div>
where $$\kappa $$ is the (real) wave number and $$\omega$$ is the (complex) frequency.  The
relation between the wave number and the frequency, $$\omega= \omega(\kappa )$$,
is called a **dispersion relation**. 

#### The von Neumann approach to stability analysis ####
A numerical algorithm for solving a PDE is said to be **stable** if the solution
at a given time does not blow up as the time step decreases, i.e., approximation errors
do not increase with time.

The von Neumann technique, proposed by John von Neumann in the 1950s expanding upon earlier
work by Crank and Nicolson, is the most popular approach for finding stability conditions
for **linear PDEs with periodic boundary conditions**.  No general approach is for stability
analysis is yet known for nonlinear PDEs.  However, for many practical nonlinear PDEs, the
linearized form is amenable to the same type of analysis.

The von Neumann approach uses the observation that linear PDEs admit plane wave solutions
to find a bound on numerical errors produced by explicit one-step finite differences.
Let us examine the technique in a bit more detail for the first-order wave equation
with forward differences.  We have

<div>
$$
  \frac{u(x, t+\Delta t) - u(x,t)}{\Delta t}
  + c\,\frac{u(x+\Delta x, t) - u(x,t)}{\Delta x} = 0 
$$
</div>
We can write the above as
<div>
$$
  u_{i,j+1} - u_{i,j} + \frac{c \Delta t}{\Delta x} (u_{i+1,j} - u_{i,j}) = 0
$$
</div>
or,
<div>
$$
  u_{i,j+1} =  (1 + C) u_{i,j} - C u_{i+1,j} \quad \text{where} \quad
  C := \frac{c \Delta t}{\Delta x} 
$$
</div>
Since the PDE admits plane wave solutions, let us examine what happens when we insert solutions
of that form into the discretized equation, i.e., 
<div>
$$
  u(x_i, t_j) = e^{i(\kappa x_i - \omega t_j)} 
$$
</div>
We then get
<div>
$$
  u_{i,j+1} = (1 + C) e^{i(\kappa x_i - \omega t_j)}
              - Ce^{i(\kappa x_{i+1} - \omega t_j)} 
$$
</div>
or
<div>
$$
  \begin{align}
  u_{i,j+1} & = (1 + C) e^{i(\kappa x_i - \omega t_j)}
              = Ce^{i(\kappa x_i - \omega t_j)} e^{i\kappa\Delta x} \\
            & = (1 + C - Ce^{i\kappa\Delta x}) u_{i,j} 
  \end{align}
$$
</div>
Therefore the expected numerical solution at $$x_i$$ is of the form
<div>
$$
  u_{j+1} = \mathcal{G}[u_{j}] 
$$
</div>
However, due to cumulative numerical rounding errors, we get a sequence of solutions
<div>
$$
  u_0 + \epsilon_0, u_1 + \epsilon_1, u_2 + \epsilon_2, \dots
$$
</div>
and the actual numerical solution that we get is
<div>
$$
  u_{j+1} + \epsilon_{j+1} = \mathcal{G}[u_j + \epsilon_j] 
$$
</div>
If we assume that a Taylor series expansion of the function $$\mathcal{G}$$ is possible,
we can make the first-order approximation
<div>
$$
  u_{j+1} + \epsilon_{j+1} = \mathcal{G}[u_j] + \epsilon_j\,\frac{\partial\mathcal{G}}{\partial u_j}
    = u_{j+1} + G\,\epsilon^j \quad \implies \epsilon_{j+1} = G\,\epsilon_j 
$$
</div>
The quantity $$G$$ is called the **amplification factor**.  For the first-order wave equation
with forward differences, from the plane wave solution we have
<div>
$$
  G = \frac{\partial\mathcal{G}}{\partial u_j}
    = \frac{\partial}{\partial u_j}\left[\left(1 + C - Ce^{i\kappa\Delta x}\right) u_{j}\right]  
    = 1 + C - Ce^{i\kappa\Delta x} 
$$
</div>
Clearly, for the solution to be stable, the amplification factor must satisfy the requirement that
<div>
$$
  |G| = \Bigl|\frac{\epsilon_{j+1}}{\epsilon_j}\Bigr| \le 1
$$
</div>
{:.notice--info}
which implies that errors do not grow in an unbounded manner with time.  For the first-order wave equation,
<div>
$$
  \begin{align}
  |G|^2 & = G\bar{G} = \left[1 + C(1-\cos\kappa\Delta x) - iC\sin\kappa\Delta x\right]
                     \left[1 + C(1-\cos\kappa\Delta x) + iC\sin\kappa\Delta x\right] \\
        & = 1 + 2C(1+C)(1-\cos\kappa\Delta x) \le 1
  \end{align}
$$
</div>
Since $$\kappa\Delta x$$ is arbitrary, the above condition is satisfied only if
<div>
$$
  2C(1+C) \le 0 \quad \implies \quad -1 \le C \le 0
   \equiv -1 \le \frac{c \Delta t}{\Delta x} \le 0
$$
</div>
You can easily show that this is equivalent to the CFL condition for this particular
backward propagating one-dimensional wave propagation problem.

For the forward propagating second-order wave equation, the equivalent stability condition is,
(from both CFL analysis and the von Neumann analysis)
<div>
$$
   0 \lt \frac{c \Delta t}{\Delta x} \le 1 \quad \equiv \quad
   \Delta t \le \frac{\Delta x}{c}
$$
</div>
{:.notice--info}

#### Remarks ####
The CFL procedure is inadequate for closed form stability bounds for many PDEs even though
convergence results can be obtained.  The von Neumann approach solves some of these issues
for linear PDEs.  

We will see how these ideas are applied to discrete element methods in the next part of this
series.

