---
layout: posts
title:  "Forward vs. Backward Euler: Plane stress plasticity"
subheadline: "Biswajit Banerjee"
description: "Part 5 of the series on plane stress Drucker-Prager plasticity"
date:  2017-03-25 10:30:00
categories:
    - Mechanics
    - Plasticity
    - Algorithm
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

##### Introduction #####
One of the main points of divergence of many implementations of plastic
return algorithms is the choice of the algorithm used for numerical integration.
Typically, this choice is limited to [forward Euler vs. backward Euler](https://en.wikipedia.org/wiki/Euler_method) for simplicity.

Recall from [Part 2]({{site.url }}/mechanics/plasticity/algorithm/plane-stress-drucker-prager-return/) that the rate equations that we need to integrate are:
<div>
$$
  \dot{\boldsymbol{\sigma}} = \mathbf{C}\dot{\boldsymbol{\varepsilon}}^e
   = \mathbf{C}\left(
    \dot{\boldsymbol{\varepsilon}} - \dot{\boldsymbol{\varepsilon}}^p\right)
$$
</div>
and
<div>
$$
\dot{\boldsymbol{\varepsilon}}^p = \dot{\lambda}\,\boldsymbol{n} = 
    \dot{\lambda} \left(
     \tfrac{1}{\sqrt{2}}
       \frac{\mathbf{P}\,\boldsymbol{\sigma}}{
         \sqrt{\boldsymbol{\sigma}^T\,\mathbf{P}\,\boldsymbol{\sigma}}
       }
     - \tfrac{1}{3} \frac{dq}{dp}\,\mathbf{I}
   \right)
$$
</div>

Let us examine way of integrating these rate equations.

##### Forward difference for stress rate #####
The first equation is typically solved using a first-order 
[forward finite difference](https://en.wikipedia.org/wiki/Finite_difference) approach:
<div>
$$
  \frac{\boldsymbol{\sigma}_{n+1} - \boldsymbol{\sigma}_n}{\Delta t} = 
  \mathbf{C}\left[
  \frac{\boldsymbol{\varepsilon}_{n+1} - \boldsymbol{\varepsilon}_n}{\Delta t}
  - \frac{\boldsymbol{\varepsilon}_{n+1}^p - \boldsymbol{\varepsilon}_n^p}{\Delta t}\right]
$$
</div>
To simplify things, a trial stress state is defined using the rate equation
<div>
$$
  \dot{\boldsymbol{\sigma}}^{\text{trial}} = \mathbf{C}\dot{\boldsymbol{\varepsilon}}
$$
</div>
and this equation is also integrated using a first-order forward difference:
<div>
$$
  \frac{\boldsymbol{\sigma}_{n+1}^{\text{trial}} - \boldsymbol{\sigma}_n^{\text{trial}}}{\Delta t} = 
  \mathbf{C}\left[
  \frac{\boldsymbol{\varepsilon}_{n+1} - \boldsymbol{\varepsilon}_n}{\Delta t}
  \right]
$$
</div>
We can then write
<div>
$$
  \frac{\boldsymbol{\sigma}_{n+1} - \boldsymbol{\sigma}_n}{\Delta t} = 
  \frac{\boldsymbol{\sigma}_{n+1}^{\text{trial}} - \boldsymbol{\sigma}_n^{\text{trial}}}{\Delta t} - 
  \mathbf{C}\left[
   \frac{\boldsymbol{\varepsilon}_{n+1}^p - \boldsymbol{\varepsilon}_n^p}{\Delta t}\right]
$$
</div>
Setting $$\boldsymbol{\sigma}_n^{\text{trial}} = \boldsymbol{\sigma}_n$$, we get the plastic return expression
<div>
$$
  \boldsymbol{\sigma}_{n+1} = 
  \boldsymbol{\sigma}_{n+1}^{\text{trial}} - 
  \mathbf{C}\left[
   \boldsymbol{\varepsilon}_{n+1}^p - \boldsymbol{\varepsilon}_n^p\right]
$$
</div>
Note that this approach is valid only if $$\mathbf{C}$$ is constant.
For a varying stiffness matrix, we have to pick either a forward Euler or a backward Euler approach to retain first-order accuracy in $$\Delta t$$.
{:.notice--info}

##### Forward and Backward Euler for stress rate #####
In many constitutive models, we get a stress-rate equation of the form
<div>
$$
  \dot{\boldsymbol{\sigma}} = \mathbf{C}(\boldsymbol{\sigma})\,\dot{\boldsymbol{\varepsilon}}^e
   = \mathbf{C}(\boldsymbol{\sigma})\left(
    \dot{\boldsymbol{\varepsilon}} - \dot{\boldsymbol{\varepsilon}}^p\right)
$$
</div>
In those situations we can use two variations on the first-order forward
difference scheme: forward Euler and backward Euler.  The forward Euler
scheme gives us the discrete form
<div>
$$
  \boldsymbol{\sigma}_{n+1} = 
  \boldsymbol{\sigma}_{n+1}^{\text{trial}} - 
  \mathbf{C}_n\left[
   \boldsymbol{\varepsilon}_{n+1}^p - \boldsymbol{\varepsilon}_n^p\right]
$$
</div>
while the backward Euler scheme leads to
<div>
$$
  \boldsymbol{\sigma}_{n+1} = 
  \boldsymbol{\sigma}_{n+1}^{\text{trial}} - 
  \mathbf{C}_{n+1}\left[
   \boldsymbol{\varepsilon}_{n+1}^p - \boldsymbol{\varepsilon}_n^p\right]
$$
</div>
Note that for problems where $$\mathbf{C}$$ is a function of the stress/deformation state, tangent modulus calculations needed by Backward Euler methods can easily become intractable and prone to bugs.  Forward Euler integration is therefore preferred for complex constitutive models, particular when there is elastic-plastic coupling.
{:.notice--warning}

##### Forward and Backward Euler for flow rule #####
For the flow rule, if we use a forward Euler scheme, we have
<div>
$$
  \frac{\boldsymbol{\varepsilon}_{n+1}^p - \boldsymbol{\varepsilon}_{n}^p}{\Delta t} = \frac{\lambda_{n+1} - \lambda_n}{\Delta t} \, \boldsymbol{n}_{n} 
$$
</div>
or
<div>
$$
  \boldsymbol{\varepsilon}_{n+1}^p = \boldsymbol{\varepsilon}_{n}^p
   + \Delta\lambda \, \boldsymbol{n}_{n} 
$$
</div>
The equivalent backward Euler update leads to
<div>
$$
  \boldsymbol{\varepsilon}_{n+1}^p = \boldsymbol{\varepsilon}_{n}^p
   + \Delta\lambda \, \boldsymbol{n}_{n+1} 
$$
</div>

##### Forward/Backward Euler stress updates #####
Using the expressions from the previous two sections, we see that for
forward Euler,
<div>
$$
  \boldsymbol{\sigma}_{n+1} = 
  \boldsymbol{\sigma}_{n+1}^{\text{trial}} - 
  \Delta\lambda\mathbf{C}_n \boldsymbol{n}_n
$$
</div>
where $$\Delta \lambda$$ can be solved using (see [Part 4]({{site.url }}/mechanics/plasticity/algorithm/plane-stress-drucker-prager-return-part-4/))
<div>
$$
  (\boldsymbol{\sigma}^\text{trial}_{n+1})^T\mathbf{P}\boldsymbol{\sigma}^\text{trial}_{n+1}
  - 2\Delta\lambda (\boldsymbol{\sigma}^\text{trial}_{n+1})^T\mathbf{P}\mathbf{C}_n\boldsymbol{n}_n
  + (\Delta\lambda)^2 (\mathbf{C}\boldsymbol{n}_n)^T\mathbf{P}\mathbf{C}_n\boldsymbol{n}_n \\
   = 2q^2(\text{tr}[\boldsymbol{\sigma}_{n+1}^\text{trial} - \Delta\lambda\mathbf{C}_n\boldsymbol{n}_n])
$$
</div>
On the other hand, for backward Euler,
<div>
$$
  \boldsymbol{\sigma}_{n+1} = 
  \boldsymbol{\sigma}_{n+1}^{\text{trial}} - 
  \Delta\lambda\mathbf{C}_{n+1} \boldsymbol{n}_{n+1}
$$
</div>
In the particular case of plane stress Drucker-Prager plasticity with
a constant elastic stiffness, we have the backward Euler scheme
<div>
$$
  \boldsymbol{\sigma}_{n+1} = 
  \boldsymbol{\sigma}_{n+1}^{\text{trial}} - 
  \Delta\lambda\mathbf{C}
    \left(
     \tfrac{1}{\sqrt{2}}
       \frac{\mathbf{P}\,\boldsymbol{\sigma}_{n+1}}{
         \sqrt{\boldsymbol{\sigma}_{n+1}^T\,\mathbf{P}\,\boldsymbol{\sigma}_{n+1}}
       }
     - \tfrac{1}{3} \frac{dq}{dp}(p_{n+1})\,\mathbf{I}
   \right)
$$
</div>
Note that in this case there is no straightforward way to cast this equation
such that we can compute $$\boldsymbol{\sigma}_{n+1}$$ in terms of a factor that scales $$\boldsymbol{\sigma}_{n+1}^{\text{trial}}$$ as is usually done for von Mises plasticity.
{:.notice--info}
Note also that for the backward Euler case, there is no straightforward
quadratic equation in $$\Delta\lambda$$ that can be solved directly.  Instead,
we need a closest point projection algorithm (or some other similar
algorithm) to find $$\Delta\lambda$$ and $$\boldsymbol{\sigma}_{n+1}$$
simultaneously.  This needs the solution of a system of equations at each
step of a Newton iteration method.
{:.notice--warning}

##### Does the choice of Forward/Backward Euler matter? #####
The forward Euler approach is easier to implement compared to the
backward Euler approach.  This is because the backward Euler approach
leads to an "implicit" set of equations that have to cast into matrix
form solved (usually iteratively) for $$\boldsymbol{\sigma}_{n+1}$$ 
and $$\Delta\lambda$$.

A question that naturally arises at this stage is whether the extra effort
needed for a backward Euler approach is justified.
{:.notice}

###### Accuracy ######
Since both approaches are first-order accurate in $$\Delta h$$, and the
backward Euler approach is typically chosen for stability in the face of
larger values of $$\Delta h$$, a stable forward Euler approach produces as
accurate solutions as the backward Euler approach for the same timestep size.

###### Stability ######
The backward Euler approach is unconditionally stable, while the stability
of the forward Euler method is limited by step size (particularly for stiff
ODE systems).  In general, the stiffness properties of complicated
constitutive equations in plasticity are non-trivial to determine as if
the safe step size.  Therefore, it is safer to use backward Euler.  However,
there is an accuracy penalty is large step sizes and backward Euler are used
and the predicted return point can deviate significantly from the exact solution.

###### Optimization ######
From an optimization point of view, the equations produced by the backward
Euler approach can be interpreted as a convex optimization problem and solved
and examined using techniques from that field.

The interpretation as an optimization problem and the robustness
of backward Euler have made it popular in commercial mechanics codes.  But
one should always keep in mind the fact that accuracy is often lost in the
process (because larger timesteps are taken). 
{:.notice--warning}

A large amount of effort is typically spent deriving consistent tangents
in implicit numerical algorithms for faster Newton-type solves.
Often, the result of these efforts is
rapid and stable convergence to inaccurate solutions.  The fact that
inaccurate solutions can have strong effects on the predictions of plasticity
models is typically ignored.
{:.notice--info}

#### Remarks ####
We mentioned the "closest-point return" approach in passing in this article.  In the next part of this series we will examine that approach in more detail.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

