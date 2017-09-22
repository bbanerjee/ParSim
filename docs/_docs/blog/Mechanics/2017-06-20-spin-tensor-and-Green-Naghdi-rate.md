---
layout: docs
title: "The difference between the spin and angular velocity tensors"
description:  "The spin tensor and the Green-Naghdi rate"
subheadline: "Biswajit Banerjee"
date:  2017-06-20 10:30:00
categories:
    - Mechanics
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

#### Introduction ####
Several years ago I took created a Wikpedia article on [objective stress rates](https://en.wikipedia.org/wiki/Objective_stress_rate) based on [lecture notes](https://en.wikiversity.org/wiki/Nonlinear_finite_elements/Objective_stress_rates) in nonlinear finite elements for a class that I taught at the University of Utah.  In this article I will attempt to clarify one recurring question on the subject: what is the difference between the spin tensor ($$\boldsymbol{w}$$) and the angular velocity tensor ($$\boldsymbol{\Omega}$$)?

#### The spin tensor ####
Following standard continuum mechanics notation, let us denote a deformation by
<div>
$$
  \mathbf{x} = \boldsymbol{\chi}(\mathbf{X}, t)
$$
</div>
where $$\mathbf{x}$$ is the position at time $$t$$ and $$\mathbf{X}$$ is the reference position of a
material point in a body.

The motion is a **rigid body translation or rotation** (or a combination of rigid body motions) if the
distance, $$D(t)$$, between two points in the body does not change with time, i.e.,
<div>
$$
  \dot{D} = \frac{d}{dt} D(t) = \frac{d}{dt} \lVert \mathbf{x}_1(t) - \mathbf{x}_2(t) \rVert = 
  \frac{\partial}{\partial t} \lVert \boldsymbol{\chi}(\mathbf{X}_1, t) - \boldsymbol{\chi}(\mathbf{X}_2, t) \rVert = 0
$$
</div>
for all points with initial positions $$\mathbf{X}_1$$ and $$\mathbf{X}_2$$.

Taking the derivative, we have
<div>
$$
  \frac{\partial}{\partial t} \lVert \boldsymbol{\chi}(\mathbf{X}_1, t) - \boldsymbol{\chi}(\mathbf{X}_2, t) \rVert
  = \frac{\boldsymbol{\chi}(\mathbf{X}_1,t) - \boldsymbol{\chi}(\mathbf{X}_2,t)}{\lVert\boldsymbol{\chi}(\mathbf{X}_1,t) - \boldsymbol{\chi}(\mathbf{X}_2,t)\rVert}
    \cdot\left[ \frac{\partial }{\partial t}\boldsymbol{\chi}(\mathbf{X}_1, t) -
           \frac{\partial }{\partial t}\boldsymbol{\chi}(\mathbf{X}_2, t) \right]
$$
</div>
or
<div>
$$
  D\dot{D} = (\mathbf{x}_1 - \mathbf{x}_2)\cdot[\mathbf{v}(\mathbf{x}_1,t) - \mathbf{v}(\mathbf{x}_2,t)]
$$
</div>
where $$\mathbf{v}$$ is the velocity.  Therefore, for a rigid body motion,
<div>
$$
  (\mathbf{x}_1 - \mathbf{x}_2)\cdot[\mathbf{v}(\mathbf{x}_1,t) - \mathbf{v}(\mathbf{x}_2,t)] = 0
$$
</div>
{:.notice}

Let us now take the gradient of this relation with respect to $$\mathbf{x}_1$$.  Then
<div>
$$
  \boldsymbol{\nabla}_{x_1} \left[
  (\mathbf{x}_1 - \mathbf{x}_2)\cdot[\mathbf{v}(\mathbf{x}_1,t) - \mathbf{v}(\mathbf{x}_2,t)]\right]
  = \mathbf{v}(\mathbf{x}_1,t) - \mathbf{v}(\mathbf{x}_2,t) +
    (\mathbf{x}_1 - \mathbf{x}_2) \cdot \boldsymbol{\nabla}_{x_1}\mathbf{v}(\mathbf{x}_1,t)
  = 0
$$
</div>
Next, if we compute the gradient with respect to $$\mathbf{x}_2$$, we have
<div>
$$
  -\boldsymbol{\nabla}_{x_2}\mathbf{v}(\mathbf{x}_2,t) - 
     [\boldsymbol{\nabla}_{x_1}\mathbf{v}(\mathbf{x}_1,t)]^T = 0
$$
</div>
If we take the limit as $$\mathbf{x}_2 \rightarrow \mathbf{x}_1 = \mathbf{x}$$, we find that
the velocity gradient, $$\boldsymbol{l}$$,  is
<div>
$$
  \boldsymbol{l} := \boldsymbol{\nabla}\mathbf{v} = -[\boldsymbol{\nabla}\mathbf{v}]^T = -\boldsymbol{l}^T \,.
$$
</div>
{:.notice}
The velocity gradient for deformable bodies is typically decomposed into a symmetric part, $$\boldsymbol{d}$$ and a skew-symmetric part $$\boldsymbol{w}$$. Clearly, from the above relation, the symmetric part is
zero and we are left with
<div>
$$
  \boldsymbol{l} = \boldsymbol{w} \,.
$$
</div>
{:.notice--info}
Therefore, the velocity gradient during a rigid body motion is called the **spin** tensor.
{:.notice--warning}

#### The Green-Naghdi objective rate ####
The Lie derivative of the Cauchy stress is given by
<div>
$$
  \overset{\circ}{\boldsymbol{\sigma}} = J^{-1}~\boldsymbol{F}\cdot
       \left[\cfrac{\partial}{\partial t}\left(J~\boldsymbol{F}^{-1}\cdot\boldsymbol{\sigma}\cdot\boldsymbol{F}^{-T}\right)\right]
       \cdot\boldsymbol{F}^T 
$$
</div>
where $$\boldsymbol{F}$$ is the deformation gradient and $$J = \det\boldsymbol{F}$$.
From the polar decomposition theorem we have
<div>
$$
  \boldsymbol{F} = \boldsymbol{R}\cdot\boldsymbol{U}
$$
</div>
where $$\boldsymbol{R}$$ is the orthogonal rotation tensor ($$\boldsymbol{R}^{-1} = \boldsymbol{R}^T$$)
and $$\boldsymbol{U}$$ is the symmetric, positive definite, right stretch.

If we assume that $$\boldsymbol{U} = \boldsymbol{\mathit{1}}$$ we get
$$\boldsymbol{F} = \boldsymbol{R}$$.  Also, since there is no stretch, $$J = 1$$.  Note that this
doesn't mean that there is not stretch in the actual body; this simplification is just
for the purposes of defining an objective stress rate.  Therefore,
<div>
$$
  \overset{\circ}{\boldsymbol{\sigma}} = \boldsymbol{R}\cdot
       \left[\cfrac{\partial}{\partial t}\left(\boldsymbol{R}^{-1}\cdot\boldsymbol{\sigma}\cdot\boldsymbol{R}^{-T}\right)\right]
       \cdot\boldsymbol{R}^T 
    = \boldsymbol{R}\cdot\left[\cfrac{\partial }{\partial t}\left(\boldsymbol{R}^T\cdot\boldsymbol{\sigma}\cdot\boldsymbol{R}\right)\right]
       \cdot\boldsymbol{R}^T 
$$
</div>
If we define the **angular velocity tensor**, $$\boldsymbol{\Omega}$$, as
<div>
$$
  \boldsymbol{\Omega} := \dot{\boldsymbol{R}}\cdot\boldsymbol{R}^T
$$
</div>
{:.notice--info}
the Lie derivative of the Cauchy stress reduces to the **Green–Naghdi** rate:
<div>
$$
  \overset{\square}{\boldsymbol{\sigma}} = \dot{\boldsymbol{\sigma}} + \boldsymbol{\sigma}\cdot\boldsymbol{\Omega}
    - \boldsymbol{\Omega}\cdot\boldsymbol{\sigma} \,.
$$
</div>
{:.notice}

#### Relation between spin and angular velocity ####
We know that
<div>
$$
  \boldsymbol{l} = \dot{\boldsymbol{F}}\cdot\boldsymbol{F}^{-1} ~,~~
  \boldsymbol{w} = \tfrac{1}{2}(\boldsymbol{l} - \boldsymbol{l}^T) \,.
$$
</div>
Therefore
<div>
$$
  \boldsymbol{w} = \tfrac{1}{2}(\dot{\boldsymbol{F}}\cdot\boldsymbol{F}^{-1} - \boldsymbol{F}^{-T}\cdot\dot{\boldsymbol{F}}^T)
$$
</div>
Using
<div>
$$
  \boldsymbol{F} = \boldsymbol{R}\cdot\boldsymbol{U} 
  \quad \implies \quad \dot{\boldsymbol{F}} = \dot{\boldsymbol{R}}\cdot\boldsymbol{U} + \boldsymbol{R}\cdot\dot{\boldsymbol{U}}
$$
</div>
we have
<div>
$$
  \dot{\boldsymbol{F}}\cdot\boldsymbol{F}^{-1} =  (\dot{\boldsymbol{R}}\cdot\boldsymbol{U} + \boldsymbol{R}\cdot\dot{\boldsymbol{U}})\cdot
     (\boldsymbol{U}^{-1}\cdot\boldsymbol{R}^T)
     = \dot{\boldsymbol{R}}\cdot\boldsymbol{R}^T + \boldsymbol{R}\cdot\dot{\boldsymbol{U}}\cdot\boldsymbol{U}^{-1}\cdot\boldsymbol{R}^T
$$
</div>
and
<div>
$$
  \boldsymbol{F}^{-T}\cdot\dot{\boldsymbol{F}}^T =  (\boldsymbol{R}\cdot\boldsymbol{U}^{-1})\cdot
      (\boldsymbol{U}\cdot\dot{\boldsymbol{R}}^T + \dot{\boldsymbol{U}}\cdot\boldsymbol{R}^T)
    = \boldsymbol{R}\cdot\dot{\boldsymbol{R}}^T + \boldsymbol{R}\cdot\boldsymbol{U}^{-1}\cdot\dot{\boldsymbol{U}}\cdot\boldsymbol{R}^T \,.
$$
</div>
Combining the two relations above in the expression for $$\boldsymbol{w}$$ gives us
<div>
$$
  \boldsymbol{w} = \frac{1}{2}~(\dot{\boldsymbol{R}}\cdot\boldsymbol{R}^T + \boldsymbol{R}\cdot\dot{\boldsymbol{U}}\cdot\boldsymbol{U}^{-1}\cdot\boldsymbol{R}^T
         - \boldsymbol{R}\cdot\dot{\boldsymbol{R}}^T - \boldsymbol{R}\cdot\boldsymbol{U}^{-1}\cdot\dot{\boldsymbol{U}}\cdot\boldsymbol{R}^T) \,.
$$
</div>
We can simplify the above relation by noting that
<div>
$$
  \boldsymbol{R}\cdot\boldsymbol{R}^T = \boldsymbol{\mathit{1}} 
  \quad \implies \quad \dot{\boldsymbol{R}}\cdot\boldsymbol{R}^T + \boldsymbol{R}\cdot\dot{\boldsymbol{R}}^T = \boldsymbol{\mathit{0}}\,.
$$
</div>
Therefore
<div>
$$
  \boldsymbol{w} = \dot{\boldsymbol{R}}\cdot\boldsymbol{R}^T + \frac{1}{2}~\boldsymbol{R}\cdot(\dot{\boldsymbol{U}}\cdot\boldsymbol{U}^{-1} - 
             \boldsymbol{U}^{-1}\cdot\dot{\boldsymbol{U}})\cdot\boldsymbol{R}^T
$$
</div>
For pure **rigid body motion**, the stretch is identity and its rate of change is zero.  So we have
<div>
$$
  \boldsymbol{w} = \dot{\boldsymbol{R}}\cdot\boldsymbol{R}^T = \boldsymbol{\Omega}
$$
</div>
{:.notice--info}
That implies that the spin tensor and the angular velocity tensor are identical to the
velocity gradient for rigid body motions.
{:.notice--warning}

#### Remarks ####
In many practical applications the spin and the angular velocity tensors are assumed to be identical.
One has to be careful to make sure that the assumptions made during that identification are valid
for the application under consideration.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

