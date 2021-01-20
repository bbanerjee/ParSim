---
layout: posts
title:  "Nonlinear programming and closest point return plasticity"
subheadline: "Biswajit Banerjee"
description: "Part 6 of the series on plasticity return algorithms"
date:  2017-03-27 10:30:00
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
In [Part 5]({{site.baseurl }}/mechanics/plasticity/algorithm/plane-stress-drucker-prager-return-part-5/) I
briefly hinted at the closest-point return algorithm.  The ideas behind this were made rigorous in the
mid-to-late 1980s by a group of researchers influenced by developments in convex optimization. Since
then there has been a statis in the development of return algorithms for phenomenological plasticity.
However, the recent surge in improvements in machine learning has the potential of leading to another
leap forward in our understanding and implementation of nonlinear, history-dependent, constitutive models.

Let us explore the basic ideas behind closest-point algorithms.

##### Background #####
In nonlinear optimization, the method of [Lagrange multipliers](https://en.wikipedia.org/wiki/Lagrange_multiplier) has been used since the mid 1800s to solve minimization problems with *equality* constraints.  In
1950, this approach was generalized by [Kuhn and Tucker](https://projecteuclid.org/euclid.bsmsp/1200500249) to allow for *inequality* constraints.  Later it was discovered that W. Karush from the University of Chicago had reached the same conclusions in his MSc thesis from 1939.

###### Primal form ######
The primal form of the optimization problem is
<div>
$$
  \begin{align}
    & \text{minimize}   & & f(\mathbf{x}) \\
    & \text{subject to} & & g_i(\mathbf{x}) \le 0, \quad i = 1, \dots, m \\
    &                   & & h_j(\mathbf{x}) = 0, \quad j = 1, \dots, p 
  \end{align}
$$
</div>
Note that there is no convexity requirement for this problem.

###### The Lagrangian ######
The Lagrangian ($$\mathcal{L}$$) associated with the primal form is just the weighted sum of
the objective function $$f_0$$ and the constraint functions $$g_i$$ and $$h_j$$.
Thus
<div>
$$
  \mathcal{L}(\mathbf{x}, \boldsymbol{\lambda}, \boldsymbol{\nu})
    = f(\mathbf{x}) + \boldsymbol{\lambda}\cdot\mathbf{g}(\mathbf{x}) +
      \boldsymbol{\nu}\cdot\mathbf{h}(\mathbf{x})
$$
</div>
where
<div>
$$
  \boldsymbol{\lambda} = \begin{bmatrix} \lambda_1 \\ \lambda_2 \\ \vdots \\ \lambda_m \end{bmatrix} ~,~~
  \mathbf{g} = \begin{bmatrix} g_1 \\ g_2 \\ \vdots \\ g_m \end{bmatrix} ~,~~
  \boldsymbol{\nu} = \begin{bmatrix} \nu_1 \\ \nu_2 \\ \vdots \\ \nu_p \end{bmatrix} ~,~~
  \mathbf{h} = \begin{bmatrix} h_1 \\ h_2 \\ \vdots \\ h_p \end{bmatrix} \,.
$$
</div>
The vectors $$\boldsymbol{\lambda}$$ and $$\boldsymbol{\nu}$$ are called *Lagrange multiplier vectors* or,
more frequently, the *dual variables* of the primal problem.

###### Dual function ######
The dual function ($$\mathcal{F}(\boldsymbol{\lambda},\boldsymbol{\nu})$$)
to the primal problem is defined as
<div>
$$
  \mathcal{F}(\boldsymbol{\lambda},\boldsymbol{\nu}) = \inf_{\mathbf{x}}
  \mathcal{L}(\mathbf{x}, \boldsymbol{\lambda}, \boldsymbol{\nu})
    = \inf_{\mathbf{x}} \left[f(\mathbf{x}) + \boldsymbol{\lambda}\cdot\mathbf{g}(\mathbf{x}) +
      \boldsymbol{\nu}\cdot\mathbf{h}(\mathbf{x})\right]
$$
</div>
Note that the dual function is the minimum of a family of affine functions (linear + a constant term)
in $$(\boldsymbol{\lambda},\boldsymbol{\nu})$$.  This makes the dual problem concave.  Note also that
since the dual function is affine, it is bounded from below by $$-\infty$$ when the value of $$\mathbf{x}$$
is unbounded.

Simplified forms for $$\mathcal{F}$$ can be found for many problems, including problems that
can be expressed as quadratic forms.

###### Dual form ######
Since the dual function is the largest lower bound on the Lagrangian, the *Lagrange dual form* of the
primal minimization can be expressed as
<div>
$$
  \begin{align}
    & \text{maximize}   & & \mathcal{F}(\boldsymbol{\lambda},\boldsymbol{\nu}) \\
    & \text{subject to} & & \boldsymbol{\lambda} \succeq \mathbf{0}
  \end{align}
$$
</div>
We don't have any constraint on $$\boldsymbol{\nu}$$ because $$\mathbf{h}(\mathbf{x}) = \mathbf{0}$$.

###### Karush-Kuhn-Tucker optimality conditions ######
Let $$\mathbf{x}^\star$$ be the optimal solution for the primal problem and
let $$(\boldsymbol{\lambda}^\star, \boldsymbol{\nu}^\star)$$ be the optimal solution
of the dual problem.  When these two solutions lead to a zero duality gap, i.e.,
<div>
$$
  f(\mathbf{x}^\star) = \mathcal{F}(\boldsymbol{\lambda}^\star, \boldsymbol{\nu}^\star)
$$
</div>
the Lagrangian at that optimal point is
<div>
$$
  \mathcal{L}(\mathbf{x}^\star, \boldsymbol{\lambda}^\star, \boldsymbol{\nu}^\star)
    = f(\mathbf{x}^\star) + \boldsymbol{\lambda}^\star\cdot\mathbf{g}(\mathbf{x}^\star) +
      \boldsymbol{\nu}^\star\cdot\mathbf{h}(\mathbf{x}^\star)
$$
</div>
Also, since $$\boldsymbol{\lambda}^\star \ge 0$$ and $$\mathbf{h} = 0$$, 
<div>
$$
  f(\mathbf{x}^\star) = \mathcal{F}(\boldsymbol{\lambda}^\star, \boldsymbol{\nu}^\star)
  = \inf_{\mathbf{x}} \mathcal{L}(\mathbf{x}, \boldsymbol{\lambda}^\star, \boldsymbol{\nu}^\star)
  \le \mathcal{L}(\mathbf{x}^\star, \boldsymbol{\lambda}^\star, \boldsymbol{\nu}^\star)
  \le f(\mathbf{x}^\star)
$$
</div>
The only way for the above to be true is when 
<div>
$$
  \boldsymbol{\lambda}^\star\cdot\mathbf{g}(\mathbf{x}^\star) = 0  \quad \leftrightarrow \quad
  \lambda^\star_i g_i(\mathbf{x}^\star) = 0 \,.
$$
</div>
Also, since $$\mathbf{x}^\star$$ minimizes the Lagrangian, its gradient is zero at that point:
<div>
$$
  \frac{\partial}{\partial \mathbf{x}}
    \mathcal{L}(\mathbf{x}^\star, \boldsymbol{\lambda}^\star, \boldsymbol{\nu}^\star)
    = \mathbf{0} 
    = \frac{\partial f(\mathbf{x}^\star)}{\partial\mathbf{x}} +
      \boldsymbol{\lambda}^\star\cdot\frac{\partial\mathbf{g}(\mathbf{x}^\star)}{\partial\mathbf{x}} +
      \boldsymbol{\nu}^\star\cdot\frac{\partial\mathbf{h}(\mathbf{x}^\star)}{\partial\mathbf{x}}
$$
</div>
These results, along with the original constraints of the primal and dual problems, are collected
together into the *Karush-Kuhn-Tucker optimality conditions*:
<div>
$$
  \begin{align}
    & g_i(\mathbf{x}^\star) \le 0 & & 
     h_j(\mathbf{x}^\star) = 0  \\
    & \lambda_i^\star \ge 0  & & 
    \lambda^\star_i g_i(\mathbf{x}^\star) = 0 \\
    & \frac{\partial f(\mathbf{x}^\star)}{\partial\mathbf{x}} +
      \boldsymbol{\lambda}^\star\cdot\frac{\partial\mathbf{g}(\mathbf{x}^\star)}{\partial\mathbf{x}} +
      \boldsymbol{\nu}^\star\cdot\frac{\partial\mathbf{h}(\mathbf{x}^\star)}{\partial\mathbf{x}} = \mathbf{0}
  \end{align}
$$
</div>
{:.notice}

##### Similarity with plasticity #####
The plastic loading-unloading conditions are similar to the Karush-Kush-Tucker optimality conditions
in that we have
<div>
$$
  \begin{align}
    g(\boldsymbol{\sigma}) \le 0 ~,~~
    \dot{\lambda} \ge 0 ~,~~ \dot{\lambda} g(\boldsymbol{\sigma}) = 0
  \end{align}
$$
</div>
where $$g(\boldsymbol{\sigma})$$ is the yield surface constraining the values of $$\boldsymbol{\sigma}$$.
We may also interpret the flow rule as the last Karush-Kuhn-Tucker condition:
<div>
$$
  -\dot{\boldsymbol{\varepsilon}}^p + \dot{\lambda}\frac{\partial g}{\partial \boldsymbol{\sigma}} = 0 
  \quad \text{where} \quad
  -\dot{\boldsymbol{\varepsilon}}^p =: \frac{\partial f}{\partial \boldsymbol{\sigma}}
$$
</div>
and $$f(\boldsymbol{\sigma})$$ is the quantity that is minimized in the primal problem.  We can
interpret $$f$$ as the negative of the maximum plastic dissipation, i.e.,
$$
  f(\boldsymbol{\sigma}) = -\boldsymbol{\sigma}:\dot{\boldsymbol{\varepsilon}}^p \,.
$$

If we use a first-order update approach, the discretized equations for perfect plasticity are
(see [Part 5]({{site.baseurl }}/mechanics/plasticity/algorithm/plane-stress-drucker-prager-return-part-5/))
<div>
$$
  \begin{align}
    \boldsymbol{\sigma}_{n+1} & = \mathbf{C}:(\boldsymbol{\varepsilon}_{n+1} - \boldsymbol{\varepsilon}_{n+1}^p)
    = \boldsymbol{\sigma}_{n+1}^{\text{trial}} - \mathbf{C}:(\boldsymbol{\varepsilon}_{n+1}^p - \boldsymbol{\varepsilon}_{n}^p)\\
    \boldsymbol{\varepsilon}_{n+1}^p & = \boldsymbol{\varepsilon}_n^p + \Delta\lambda \left.\frac{\partial g}{\partial \boldsymbol{\sigma}}\right|_{\boldsymbol{\sigma}_{n}}
    \quad \text{or} \quad
    \boldsymbol{\varepsilon}_{n+1}^p = \boldsymbol{\varepsilon}_n^p + \Delta\lambda \left.\frac{\partial g}{\partial \boldsymbol{\sigma}}\right|_{\boldsymbol{\sigma}_{n+1}} \\
    & g(\boldsymbol{\sigma}_{n+1})  \le 0 ~,~~
    \Delta\lambda \ge 0 ~,~~ \Delta\lambda g(\boldsymbol{\sigma}_{n+1}) = 0
  \end{align}
$$
</div>
Note that if we interpret the flow rule as an optimality condition
a backward Euler update is consistent with the Karush-Kuhn-Tucker conditions
and a forward Euler update is ruled out.
{:.notice--info}

##### Closest point return #####
Let $$\boldsymbol{\sigma}^{\text{trial}}$$ be the trial stress and let
$$g(\boldsymbol{\sigma}^{\text{trial}})$$ be the value of the yield function at that state.  Let
$$\boldsymbol{\sigma}_{n+1}$$ be actual stress and let $$g(\boldsymbol{\sigma}_{n+1}) = 0$$ be the value
of the yield function at the actual stress state.

Let us assume the actual stress state on the yield surface is at the closest distance from the trial stress.
Then we can devise the primal minimization problem:
<div>
$$
  \begin{align}
    & \text{minimize}   & & f(\boldsymbol{\sigma}) = \lVert \boldsymbol{\sigma}^{\text{trial}} - \boldsymbol{\sigma}\rVert^2 \\
    & \text{subject to} & & g(\boldsymbol{\sigma}) \le 0 \\
  \end{align}
$$
</div>
where
<div>
$$
  \lVert \boldsymbol{\sigma} \rVert = \sqrt{\boldsymbol{\sigma}:\boldsymbol{\sigma}}
$$
</div>
The Lagrangian for this problem is
<div>
$$
  \mathcal{L}(\boldsymbol{\sigma},\lambda) =
  f(\boldsymbol{\sigma})+ \Delta\lambda g(\boldsymbol{\sigma}) = 
  \lVert \boldsymbol{\sigma}^{\text{trial}} - \boldsymbol{\sigma}\rVert^2 + \Delta\lambda g(\boldsymbol{\sigma})
$$
</div>
The Karush-Kuhn-Tucker conditions for this problem at the optimum value $$\boldsymbol{\sigma}_{n+1}$$ are
<div>
$$
  \begin{align}
    & g(\boldsymbol{\sigma}_{n+1}) \le 0 ~,~~
      \Delta\lambda \ge 0 ~,~~
      \Delta\lambda g(\boldsymbol{\sigma}_{n+1}) = 0 \\
    & \frac{\partial f(\boldsymbol{\sigma}_{n+1})}{\partial\boldsymbol{\sigma}} + \Delta\lambda \frac{\partial g(\boldsymbol{\sigma}_{n+1})}{\partial \boldsymbol{\sigma}} =  
      -2(\boldsymbol{\sigma}^{\text{trial}} - \boldsymbol{\sigma}_{n+1}) + \Delta\lambda \frac{\partial g(\boldsymbol{\sigma}_{n+1})}{\partial \boldsymbol{\sigma}} = \mathbf{0} 
  \end{align}
$$
</div>
From the last condition we see that the closest distance using this criterion leads to a stress value of
<div>
$$
   \boldsymbol{\sigma}_{n+1} = 
      \boldsymbol{\sigma}^{\text{trial}} -  \tfrac{1}{2} \Delta\lambda \frac{\partial g(\boldsymbol{\sigma}_{n+1})}{\partial \boldsymbol{\sigma}}
$$
</div>
But we have seen previously that the first-order stress update with backward Euler leads to
<div>
$$
  \boldsymbol{\sigma}_{n+1}  
    = \boldsymbol{\sigma}^{\text{trial}} - \Delta\lambda\mathbf{C}: \frac{\partial g(\boldsymbol{\sigma}_{n+1})}{\partial \boldsymbol{\sigma}}
$$
</div>
The similarity between the two indicates that we are on the right track, i.e., the actual stress is at
the closest distance from the trial stress to the yield surface.  But the correct closest distance is
not in the standard standard stress space, but in a space where the norm to be minimized is given by
<div>
$$
  \lVert \boldsymbol{\sigma} \rVert_{\mathbf{C}^{-1}} = \sqrt{\boldsymbol{\sigma}:\mathbf{C}^{-1}:\boldsymbol{\sigma}}
$$
</div>
{:.notice--info}

This can be verified by repeating the above exercise with the new definition of the norm.

#### Remarks ####
Most methods used for finding $$\boldsymbol{\sigma}_{n+1}$$ using the closest-point projection idea
use variations on the Newton method that require the computation of second-derivatives of the
yield function.  We will discuss a method that avoids those computations in the next part of this series.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

