---
layout: posts
title:  "Material and spatial incremental constitutive equations"
subheadline: "Biswajit Banerjee"
description: "An answer to a common question on objectivity"
date:  2017-05-31 10:30:00
categories:
    - Mechanics
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

#### The question ####
A colleague asked a question on objectivity a few days ago that had me going back to
Ray Ogden's book on nonlinear elastic deformations.  The question was on incremental
stresses and their material and spatial descriptions.

To be more specific, the question was on incremental stress-strain equations expressed
in rate form and why the instantaneous moduli for material spatial stress measures were
different.

###### Instantaneous moduli for PK-2 stress and Green strain ######
Consider the relation between the second Piola-Kirchhoff (PK-2) stress
($$\boldsymbol{S}$$) and the Lagrangian Green strain ($$\boldsymbol{E}$$):
<div>
$$
  \dot{\boldsymbol{S}} = \mathsf{C}^{SE} : \dot{\boldsymbol{E}}
$$
</div>
where $$\mathsf{C}^{SE}$$ is the first-order instantaneous modulus tensor (a rank-4 tensor).

Note that the first-order instantaneous modulus tensor is neither constant
nor independent of the state of deformation.
{:.notice--info}

Many books use an additional linearization assumption to simplify the instantaneous modulus tensor.
In particular, for an isotropic material, the linearized instantaneous modulus tensor can be
expressed as
<div>
$$
  \mathsf{C}^{SE} = \kappa \mathbf{I}\otimes\mathbf{I} +
    2\mu\left(\mathsf{I}^{(4s)} - \tfrac{1}{3}\mathbf{I}\otimes\mathbf{I}\right)
$$
</div>
where $$\kappa$$ and $$\mu$$ are material constants, $$\mathbf{I}$$ is the rank-2 identity tensor,
and $$\mathbf{I}^{(4s)}$$ is symmetrized rank-4 identity tensor.

In our [earlier article on closest-point return]({{ site.baseurl }}/mechanics/plasticity/algorithm/closest-point-return/) we saw that we can express the above modulus tensor in the terms of projection tensors:
<div>
$$
  \mathsf{C}^{SE} = 3\kappa \mathsf{P}^{\text{iso}} + 2\mu\mathsf{P}^{\text{symdev}} \quad \quad \text{(1)}
$$
</div>
{:.notice}

###### Instantaneous moduli for Kirchhoff stress ######
My colleague wanted to use the expression in equation (1) above to determine the instantaneous moduli
for a system that uses the Kirchhoff stress measure ($$\boldsymbol{\tau}$$)
<div>
$$
  \boldsymbol{\tau} = J \boldsymbol{\sigma} = \boldsymbol{F}\cdot\boldsymbol{S}\cdot\boldsymbol{F}^T
$$
</div>
where $$\boldsymbol{F}$$ is the deformation gradient, $$J = \det\boldsymbol{F}$$, and $$\boldsymbol{\sigma}$$ is
the Cauchy stress.

The time derivative of this stress measure is
<div>
$$
  \begin{align}
  \dot{\boldsymbol{\tau}} & = \frac{d}{dt}(\boldsymbol{F}\cdot\boldsymbol{S}\cdot\boldsymbol{F}^T) \\
     & = \dot{\boldsymbol{F}}\cdot\boldsymbol{S}\cdot\boldsymbol{F}^T + 
         \boldsymbol{F}\cdot\dot{\boldsymbol{S}}\cdot\boldsymbol{F}^T + 
         \boldsymbol{F}\cdot\boldsymbol{S}\cdot\dot{\boldsymbol{F}^T} \\
     & = \boldsymbol{l}\cdot\boldsymbol{F}\cdot\boldsymbol{S}\cdot\boldsymbol{F}^T + 
         \boldsymbol{F}\cdot(\mathsf{C}^{SE}:\dot{\boldsymbol{E}})\cdot\boldsymbol{F}^T + 
         \boldsymbol{F}\cdot\boldsymbol{S}\cdot\boldsymbol{F}^T\cdot\boldsymbol{l}^T
  \end{align}
$$
</div>
where we have used the rate-form expression for $$\dot{\boldsymbol{S}}$$ and the relationship
between the velocity gradient ($$\boldsymbol{l}$$) and time derivative of the deformation gradient.
Also,
<div>
$$
  \begin{align}
  \dot{\boldsymbol{E}}
    & = \tfrac{1}{2}(\dot{\boldsymbol{F}^T}\cdot\boldsymbol{F} + \boldsymbol{F}^T\cdot\dot{\boldsymbol{F}}) \\
    & = \tfrac{1}{2}(\boldsymbol{F}^T\cdot\boldsymbol{l}^T\cdot\boldsymbol{F} +
                     \boldsymbol{F}^T\cdot\boldsymbol{l}\cdot\boldsymbol{F}) \\
    & = \boldsymbol{F}^T\cdot\boldsymbol{d}\cdot\boldsymbol{F}
  \end{align}
$$
</div>
where $$\boldsymbol{d}$$ is the symmetric part of the velocity gradient tensor.
Therefore, defining the spin tensor ($$\boldsymbol{w}$$) via $$\boldsymbol{l} = \boldsymbol{d} + \boldsymbol{w}$$,
we have
<div>
  \begin{align}
  \dot{\boldsymbol{\tau}} 
     & = \boldsymbol{l}\cdot\boldsymbol{\tau} + 
         \boldsymbol{F}\cdot\mathsf{C}^{SE}:(\boldsymbol{F}^T\cdot\boldsymbol{d}\cdot\boldsymbol{F})\cdot\boldsymbol{F}^T + 
         \boldsymbol{\tau}\cdot\boldsymbol{l}^T \\
     & = \boldsymbol{d}\cdot\boldsymbol{\tau} +  \boldsymbol{w}\cdot\boldsymbol{\tau} + 
         \boldsymbol{F}\cdot\mathsf{C}^{SE}:(\boldsymbol{F}^T\cdot\boldsymbol{d}\cdot\boldsymbol{F})\cdot\boldsymbol{F}^T + 
         \boldsymbol{\tau}\cdot\boldsymbol{d}^T + \boldsymbol{\tau}\cdot\boldsymbol{w}^T 
  \end{align}
</div>
If we define the Jaumann rate of the Kirchhoff stress as
<div>
$$
  \overset{\triangle J}{\boldsymbol{\tau}} = \dot{\boldsymbol{\tau}} - \boldsymbol{w}\cdot\boldsymbol{\tau} -\boldsymbol{\tau}\cdot\boldsymbol{w}^T 
$$
</div>
and the Jaumann modulus using
<div>
$$
  \mathsf{C}^{\tau J}:\boldsymbol{d} = 
    \boldsymbol{F}\cdot\mathsf{C}^{SE}:(\boldsymbol{F}^T\cdot\boldsymbol{d}\cdot\boldsymbol{F})\cdot\boldsymbol{F}^T 
$$
</div>
{:.notice}
we have, assuming that the stress is always significantly smaller than the modulus, 
<div>
$$
  \overset{\triangle J}{\boldsymbol{\tau}} = \boldsymbol{d}\cdot\boldsymbol{\tau}  +
    \boldsymbol{\tau}\cdot\boldsymbol{d}^T  + \mathsf{C}^{\tau J}:\boldsymbol{d}
   \approx \mathsf{C}^{\tau J}:\boldsymbol{d}
$$
</div>

Now, some authors express the Jaumann modulus as
<div>
$$
  \mathsf{C}^{\tau J} = 3\kappa \mathsf{P}^{\text{iso}} + 2\mu\mathsf{P}^{\text{symdev}} \quad \quad \text{(2)}
$$
</div>
{:.notice}
**Where has the dependence on the deformation gradient gone?**  This is a common source of confusion.
{:.notice--warning}

#### The reason for the inconsistency ####
The main reason for this inconsistency is

* the use of the same symbols for two quite different sets of moduli and basis tensors, and
* ignoring the fact that a linearization operation has been performed to get the simplified modulus tensors.

In other words, most importantly, 

1. the quantities $$\kappa, \mu$$ are not necessarily the same for the two cases,
2. the basis tensors $$\mathsf{P}^{\text{iso}}$$ and $$\mathsf{P}^{\text{symdev}}$$ are not
identical for the PK-2 case and the Kirchhoff stress case.  One has been rotated and stretched relative
to the other.

A detailed discussion of these issues can be found in Odgen's book (chapter 6.1.4).  A shorter discussion
can be found in *Computational Inelasticity* by Simo and Hughes (sections 7.1.5.3 - 7.1.5.5).

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

