---
title:  "Exploring closest point return plasticity"
subheadline: "Biswajit Banerjee"
description: "Part 7 of the series on plasticity return algorithms"
date:  2017-03-28 10:30:00
categories:
    - Mechanics
    - Plasticity
    - Algorithm
excerpt_separator: <!--more-->
toc: true
toc_label: "Contents"
toc_sticky: true
---


##### Introduction #####
In [Part 6]({{site.baseurl }}/mechanics/plasticity/algorithm/primal-dual-closest-point-return/), I
explained why a backward Euler stress update and a closest point return from the trial stress
to the yield surface are closely related.  More specifically, the correct updated stress
is at the shortest distance from the trial stress to the yield surface
in a 9-dimensional space that has the Euclidean distance measure
<!--more-->
<div>
$$
  \lVert \boldsymbol{\sigma} \rVert_{\mathsf{C}^{-1}} = \sqrt{\boldsymbol{\sigma}:\mathsf{C}^{-1}:\boldsymbol{\sigma}}
$$
</div>
where $$\mathsf{C}$$ is the stiffness tensor.  We will explore some of the implications of this idea
in this article.

Note that this particular closest-point interpretation applies only for *perfect plasticity* and
only *associative* flow rules.  For hardening plasticity, the space in which the actual stress is
closest to the trial stress is different.  For non-associative plasticity, it is unclear whether any
closest-point approach can be rigorously justified.
{:.notice--warning}

##### Eigendecompositions in linear elasticity #####
The stiffness tensor for an isotropic elastic material is
<div>
$$
  \mathsf{C} = \lambda\,\mathbf{I}\otimes\mathbf{I} + 2\mu\,\mathsf{I}^s
$$
</div>
where $$\lambda,\mu$$ are the Lam&eacute; elastic constants, $$\mathbf{I}$$ is the rank-2 identity tensor, and
$$\mathsf{I}^s$$ is the symmetric rank-4 identity tensor.  The inverse of $$\mathsf{C}$$ is the
compliance tensor
<div>
$$
  \mathsf{C}^{-1} = \mathsf{S}  =
   -\frac{\lambda}{2\mu(3\lambda+2\mu)}\,\mathbf{I}\otimes\mathbf{I} + \frac{1}{2\mu}\,\mathsf{I}^s 
$$
</div>

Eigendecompositions of the stiffness and compliance tensors are defined via
<div>
$$
  \mathsf{C}:\mathbf{V} = \lambda \mathbf{V} ~,~~
  \mathsf{S}:\mathbf{V} = \frac{1}{\lambda} \mathbf{V} 
$$
</div>
where $$\lambda$$ are the eigenvalues (not to be confused with the Lam&eacute; modulus)
and $$\mathbf{V}$$ are rank-2 tensors that form the eigenbasis.  Because of the symmetries
of the stiffness matrix, there are six or less unique eigenvalues and the corresponding
eigentensors are orthogonal, i.e.,
<div>
$$
  \mathbf{V}_i : \mathbf{V}_i = 1 \quad \text{and} \quad
  \mathbf{V}_i : \mathbf{V}_j = 0  \,.
$$
</div>
The stiffness and compliance tensors may then be represented as:
<div>
$$
  \mathsf{C} = \sum_{i=1}^m \lambda_i \mathbf{V}_i\otimes\mathbf{V}_i ~,~~
  \mathsf{S} = \sum_{i=1}^m \frac{1}{\lambda_i} \mathbf{V}_i\otimes\mathbf{V}_i
$$
</div>
where $$m$$ is the number of non-zero and distinct eigenvalues.  Note also that, in this
eigenbasis, the symmetric rank-4 identity tensor is
<div>
$$
  \mathsf{I}^s = \sum_{i=1}^m \mathbf{V}_i\otimes\mathbf{V}_i \,.
$$
</div>
Eigenprojectors are defined as rank-4 tensors that have the property (for $$ i \ne j$$)
<div>
$$
  \mathsf{P}_i:\mathbf{V}_i = \mathbf{V}_i  ~,~~
  \mathsf{P}_i:\mathbf{V}_j = \mathbf{0} 
$$
</div>
If we apply the eigenprojector to the rank-4 identity tensor, we get
<div>
$$
  \mathsf{P}_k = \mathsf{P}_k : \mathsf{I}^s = 
     \sum_{i=1}^m \mathsf{P}_k: (\mathbf{V}_i\otimes\mathbf{V}_i) 
   = \sum_{i=1}^m (\mathsf{P}_k: \mathbf{V}_i)\otimes\mathbf{V}_i 
   = (\mathsf{P}_k : \mathbf{V}_k)\otimes\mathbf{V}_k
   = \mathbf{V}_k\otimes\mathbf{V}_k
$$
</div>
Therefore we may also write the eigendecomposition in terms of the eigenprojectors
<div>
$$
  \mathsf{C} = \sum_{i=1}^m \lambda_i \mathsf{P}_i ~,~~
  \mathsf{S} = \sum_{i=1}^m \frac{1}{\lambda_i} \mathsf{P}_i ~,~~
  \mathsf{I}^s = \sum_{i=1}^m \mathsf{P}_i \,.
$$
</div>
For isotropic materials, a small amount of algebra shows that there are two unique eigenvectors
which lead to the decomposition
<div>
$$
  \mathsf{C} = \lambda_1 \mathsf{P}_1 + \lambda_2 \mathsf{P}_2
  \quad \text{where} \quad
  \mathsf{S} = \frac{1}{\lambda_1} \mathsf{P}_1 + \frac{1}{\lambda_2} \mathsf{P}_2
$$
</div>
and
<div>
$$
  \mathsf{P}_1 = \tfrac{1}{3} \mathbf{I}\otimes\mathbf{I} ~,~~
  \mathsf{P}_2 = \mathsf{I}^s - \mathsf{P}_1 \,.
$$
</div>
We can now express the stiffness and compliance tensors in terms of these eigenprojections:
<div>
$$
  \begin{align}
  \mathsf{C} & = \left(\kappa-\tfrac{2}{3}\mu\right)\,\mathbf{I}\otimes\mathbf{I} +
               2\mu\,\mathsf{I}^s \\
    & = 3\kappa\left(\tfrac{1}{3}\mathbf{I}\otimes\mathbf{I}\right) +
      2\mu\,\left(\mathsf{I}^s - \tfrac{1}{3}\mathbf{I}\otimes\mathbf{I}\right)
  \end{align}
$$
</div>
where $$\kappa$$ is the bulk modulus and $$\mu$$ is the shear modulus. Also,
<div>
$$
  \begin{align}
  \mathsf{S} 
   & = \frac{1}{3}\left(\frac{1}{3\kappa} - \frac{1}{2\mu}\right)\,\mathbf{I}\otimes\mathbf{I} + \frac{1}{2\mu}\,\mathsf{I}^s \\
   & = \frac{1}{3\kappa}\left(\tfrac{1}{3}\mathbf{I}\otimes\mathbf{I}\right) +
       \frac{1}{2\mu}\left(\mathsf{I}^s - \tfrac{1}{3}\mathbf{I}\otimes\mathbf{I}\right)  
  \end{align}
$$
</div>
Therefore, we can write
<div>
$$
  \mathsf{C} = 3\kappa\,\mathsf{P}^{\text{iso}} + 2\mu\,\mathsf{P}^{\text{symdev}}
  \quad \text{and} \quad
  \mathsf{S} =  \frac{1}{3\kappa}\,\mathsf{P}^{\text{iso}} +
       \frac{1}{2\mu}\,\mathsf{P}^{\text{symdev}} 
$$
</div>
{:.notice}
where $$\mathsf{P}^{\text{iso}} = \mathsf{P}_1$$ and
$$\mathsf{P}^{\text{symdev}} = \mathsf{P}_2$$.

It is also worth noting that if
<div>
$$
  \mathsf{C}^{1/2} : \mathsf{C}^{1/2} := \mathsf{C} \quad \text{and} \quad
  \mathsf{S}^{1/2} : \mathsf{S}^{1/2} := \mathsf{S} \\
$$
</div>
then, using the property that $$\mathsf{P}_1 : \mathsf{P}_2 = \mathsf{0}$$, 
<div>
$$
  \mathsf{C}^{1/2} = \sqrt{\lambda_1} \mathsf{P}_1 + \sqrt{\lambda_2} \mathsf{P}_2
  \quad \text{where} \quad
  \mathsf{S}^{1/2} = \frac{1}{\sqrt{\lambda_1}} \mathsf{P}_1 + \frac{1}{\sqrt{\lambda_2}} \mathsf{P}_2
$$
</div>
In that case, we have
<div>
$$
  \mathsf{C}^{1/2} = \sqrt{3\kappa}\,\mathsf{P}^{\text{iso}} + \sqrt{2\mu}\,\mathsf{P}^{\text{symdev}}
  \quad \text{and} \quad
  \mathsf{S}^{1/2} =  \frac{1}{\sqrt{3\kappa}}\,\mathsf{P}^{\text{iso}} +
       \frac{1}{\sqrt{2\mu}}\,\mathsf{P}^{\text{symdev}} 
$$
</div>
{:.notice}

##### The transformed space for isotropic linear elasticity #####
Details of the transformed space for isotropic linear elasticity were worked out by M. Homel in
his 2014 PhD dissertation.  We will follow his approach in this section.

The distance measure
<div>
$$
  \lVert \boldsymbol{\sigma} \rVert_{\mathsf{S}} = \sqrt{\boldsymbol{\sigma}:\mathsf{S}:\boldsymbol{\sigma}}
$$
</div>
can be interpreted as a standard Euclidean distance measure in a transformed stress space by
observing that
<div>
$$
  \begin{align}
  \lVert \boldsymbol{\sigma} \rVert_{\mathsf{S}} & =
    \sqrt{(\boldsymbol{\sigma}:\mathsf{S}^{1/2}):(\mathsf{S}^{1/2}:\boldsymbol{\sigma})}
  \quad \text{where} \quad
  \mathsf{S}^{1/2} : \mathsf{S}^{1/2} := \mathsf{S} \\
  & = \sqrt{(\mathsf{S}^{1/2}:\boldsymbol{\sigma}):(\mathsf{S}^{1/2}:\boldsymbol{\sigma})}
  \quad \text{using the major symmetry of} \quad \mathsf{S} \\
  & = \sqrt{\boldsymbol{\sigma}^\star:\boldsymbol{\sigma}^\star}
    = \lVert \boldsymbol{\sigma}^\star \rVert
  \end{align}
$$
</div>

We would like to calculate the transformed stress tensor.

###### The Lode invariants and the Lode basis ######
The Lode basis (described by R. M. Brannon in 2009) is an alternative basis that can
be used to decompose the stress tensor.  Let us define the following deviatoric quantities:
<div>
$$
  \boldsymbol{s} = \text{dev}(\boldsymbol{\sigma})
    = \boldsymbol{\sigma} - \tfrac{1}{3}\text{tr}(\boldsymbol{\sigma})\mathbf{I}
  \quad \text{and} \quad
  \boldsymbol{t} = \text{dev}(\boldsymbol{s}\cdot\boldsymbol{s})
    = \boldsymbol{s}\cdot\boldsymbol{s} - \tfrac{1}{3}\text{tr}(\boldsymbol{s}\cdot\boldsymbol{s})\mathbf{I}
$$
</div>
The quantity $$\boldsymbol{t}$$ is also called the *Hill tensor*.

The Lode invariants of a stress tensor are
<div>
$$
  z = \tfrac{1}{\sqrt{3}}\,\text{tr}(\boldsymbol{\sigma}) ~,~~
  r = \lVert \boldsymbol{s} \rVert ~,~~
  \sin 3\theta = 3\sqrt{6}\det\left(\frac{\boldsymbol{s}}{\lVert\boldsymbol{s}\rVert}\right)
$$
</div>
These invariants are associated with an orthonormal set of unit tensors
<div>
$$
  \mathbf{E}_z = \tfrac{1}{\sqrt{3}}\,\mathbf{I} ~,~~
  \mathbf{E}_r = \frac{\boldsymbol{s}}{\lVert\boldsymbol{s}\rVert} ~,~~
  \mathbf{E}_\theta =
    \frac{\frac{\boldsymbol{t}}{\lVert\boldsymbol{t}\rVert} -
    \sin 3\theta\frac{\boldsymbol{s}}{\lVert\boldsymbol{s}\rVert}}{\cos 3\theta}
$$
</div>
The stress can be expressed in terms of the Lode basis as
<div>
$$
  \boldsymbol{\sigma} = z\,\mathbf{E}_z + r\,\mathbf{E}_r \,.
$$
</div>
{:.notice}

###### The transformed stress tensor ######
We can now compute the transformed stress tensor:
<div>
$$
  \boldsymbol{\sigma}^\star = \mathsf{S}^{1/2}:\boldsymbol{\sigma}
  = \left[\frac{1}{\sqrt{3\kappa}}\,\mathsf{P}^{\text{iso}} +
       \frac{1}{\sqrt{2\mu}}\,\mathsf{P}^{\text{symdev}} \right]:
      ( z\,\mathbf{E}_z + r\,\mathbf{E}_r) \,.
$$
</div>
We can show that
<div>
$$
  \mathsf{P}^{\text{iso}}:\mathbf{E}_z = \mathbf{E}_z ~,~~
  \mathsf{P}^{\text{iso}}:\mathbf{E}_r = \mathbf{0} ~,~~
  \mathsf{P}^{\text{symdev}}:\mathbf{E}_z = \mathbf{0} ~,~~
  \mathsf{P}^{\text{symdev}}:\mathbf{E}_r = \mathbf{E}_r 
$$
</div>
Therefore, 
<div>
$$
  \boldsymbol{\sigma}^\star 
  = \frac{z}{\sqrt{3\kappa}}\,\mathbf{E}_z +
    \frac{r}{\sqrt{2\mu}}\,\mathbf{E}_r
$$
</div>
{:.notice--info}

#### Remarks ####
So we have a straightforward way of computing stresses in the transformed space
and will use this idea in the geometrical closest point return algorithm that we
will discuss in the next part of this series.

