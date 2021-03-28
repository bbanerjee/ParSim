---
title:  "Plane stress Drucker-Prager return algorithm"
subheadline: "Biswajit Banerjee"
description: ""
date:  2017-03-04 10:30:00
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
Recently I've encountered questions on how the radial return algorithm works when applied
to plane stress plasticity.  There seems to be some confusion about the application of the
plane strain constraint.  Here's the first part of my take on the issue based on Simo's work from the mid-1980s.
<!--more-->

Warning: I have used the same symbols for tensors, matrices, and vectors (Voigt notation).  The meaning depends on the context.
{:.notice--warning}

##### Plane stress elasticity #####
In plane stress we assume that one of the eigenvalues of the stress tensor is zero.  
The standard engineering convention is to assume that the plane defined by the basis vectors
$$\mathbf{e}_1$$ and $$\mathbf{e}_2$$ contains the non-zero eigenvalues while the $$\mathbf{e}_3$$
vector indicates the direction in which the eigenvalue is zero.  Therefore,

<div>
$$
  \boldsymbol{\sigma} =
  \begin{bmatrix} \sigma_{11} & \sigma_{12} & 0 \\ \sigma_{12} & \sigma_{22} & 0 \\ 0 & 0 & 0 \end{bmatrix}
  \equiv \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{12} \end{bmatrix}
$$
</div>

Similarly, the engineering convention for the plane stress strain tensor leads to

<div>
$$
  \boldsymbol{\varepsilon} =
  \begin{bmatrix} \varepsilon_{11} & \varepsilon_{12} & 0 \\ \varepsilon_{12} & \varepsilon_{22} & 0 \\ 0 & 0 & \varepsilon_{33} \end{bmatrix}
  \equiv \begin{bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ 2\varepsilon_{12} \end{bmatrix}
$$
</div>

where we have ignored $$\varepsilon_{33}$$ because it can be expressed as a function of $$\varepsilon_{11}$$ and
$$\varepsilon_{22}$$.

For an isotropic material that is linear in its elastic range,

<div>
$$
  \boldsymbol{\sigma} = \lambda\, \text{tr}(\boldsymbol{\varepsilon})\, \mathbf{I} + 2 \mu\, \boldsymbol{\varepsilon}
$$
</div>

and we can write

<div>
$$
  \begin{bmatrix} \sigma_{11} & \sigma_{12} & 0 \\ \sigma_{12} & \sigma_{22} & 0 \\ 0 & 0 & 0 \end{bmatrix}
  = \lambda ( \varepsilon_{11} + \varepsilon_{22} + \varepsilon_{33})
    \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}
    + 2\mu \begin{bmatrix} \varepsilon_{11} & \varepsilon_{12} & 0 \\ \varepsilon_{12} & \varepsilon_{22} & 0 \\ 0 & 0 & \varepsilon_{33} \end{bmatrix}
$$
</div>

From the last row we can see that
<div>
$$
  \lambda ( \varepsilon_{11} + \varepsilon_{22} + \varepsilon_{33} ) = -2\mu\varepsilon_{33}
  \quad \implies \quad
  \varepsilon_{33} = -\frac{\lambda}{\lambda + 2\mu}( \varepsilon_{11} + \varepsilon_{22} )
$$
</div>

Therefore,

<div>
$$
\begin{align}
  \begin{bmatrix} \sigma_{11} & \sigma_{12} & 0 \\ \sigma_{12} & \sigma_{22} & 0 \\ 0 & 0 & 0 \end{bmatrix}
  & = 2\mu \begin{bmatrix} \varepsilon_{11} - \varepsilon_{33} & \varepsilon_{12} & 0 \\ \varepsilon_{12} & \varepsilon_{22} - \varepsilon_{33} & 0 \\ 0 & 0 & 0 \end{bmatrix} \\
  & = \frac{2\mu}{\lambda+2\mu} \begin{bmatrix} 2(\lambda+\mu)\varepsilon_{11} + \lambda\varepsilon_{22} & (\lambda+2\mu)\varepsilon_{12} & 0 \\ (\lambda+2\mu)\varepsilon_{12} & \lambda\varepsilon_{11}+2(\lambda+\mu)\varepsilon_{22} & 0 \\ 0 & 0 & 0 \end{bmatrix}
\end{align}
$$
</div>

In engineering notation
<div>
$$
  \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{12} \end{bmatrix} =
  \begin{bmatrix} \frac{4\mu(\lambda+\mu)}{\lambda+2\mu} & \frac{2\mu\lambda}{\lambda+2\mu} & 0 \\
                  \frac{2\mu\lambda}{\lambda+2\mu} & \frac{4\mu(\lambda+\mu)}{\lambda+2\mu} & 0 \\
                  0 & 0 & \mu
  \end{bmatrix}
  \begin{bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ 2\varepsilon_{12} \end{bmatrix}
  \quad \equiv \quad
  \boldsymbol{\sigma} = \mathbf{C}\, \boldsymbol{\varepsilon}
$$
</div>
{:.notice}

Note carefully that the stiffness matrix $$\mathbf{C}$$ is not the same as the three-dimensional
stiffness matrix for
an isotropic linear elastic material with the zero stress rows and columns removed.
{:.notice--warning}

##### Drucker-Prager plasticity #####
The Drucker-Prager yield condition can be expressed as

<div>
$$
  f(\boldsymbol{\sigma}) = \sqrt{J_2} - q(p) = 0 \quad \text{where} \quad
  J_2 = \tfrac{1}{2} \boldsymbol{s}:\boldsymbol{s} ~,~~ p = \tfrac{I_1}{3} = \tfrac{1}{3} \text{tr}(\boldsymbol{\sigma})
$$
</div>
and $$q(p)$$ is a scalar-valued function of the mean stress $$(p)$$ with
$$
  \boldsymbol{s} = \boldsymbol{\sigma} - p\, \mathbf{I} \,.
$$
The associated flow rule for the yield function $$f(\boldsymbol{\sigma})$$ is
<div>
$$
  \dot{\boldsymbol{\varepsilon}}^p = \dot{\lambda}\,\frac{\partial f}{\partial \boldsymbol{\sigma}}
$$
</div>
The derivative of $$f$$ is
<div>
$$
  \frac{\partial f}{\partial \boldsymbol{\sigma}} =
    \tfrac{1}{\sqrt{2}} \frac{\partial}{\partial \boldsymbol{\sigma}}\left[\sqrt{\boldsymbol{s}:\boldsymbol{s}}\right] - \frac{dq}{dp}\,\frac{\partial p}{\partial \boldsymbol{\sigma}}
   = \tfrac{1}{2\sqrt{2}} \frac{1}{\sqrt{\boldsymbol{s}:\boldsymbol{s}}}
     \frac{\partial}{\partial \boldsymbol{\sigma}}(\boldsymbol{s}:\boldsymbol{s}) -
     \tfrac{1}{3} \frac{dq}{dp}\,\frac{\partial \text{tr}(\boldsymbol{\sigma})}{\partial \boldsymbol{\sigma}}
$$
</div>
Now
<div>
$$
   \frac{\partial}{\partial \boldsymbol{\sigma}}(\boldsymbol{s}:\boldsymbol{s}) \equiv
   \frac{\partial}{\partial \sigma_{mn}}(s_{ij} s_{ij}) = 2 s_{ij} \frac{\partial s_{ij}}{\partial \sigma_{mn}} =
   2 s_{ij} \left[\frac{\partial \sigma_{ij}}{\partial \sigma_{mn}} - \tfrac{1}{3} \frac{\partial \sigma_{kk}}{\partial \sigma_{mn}}\,\delta_{ij}\right]
$$
</div>
or,
<div>
$$
   \frac{\partial}{\partial \boldsymbol{\sigma}}(\boldsymbol{s}:\boldsymbol{s}) \equiv
   2 s_{ij} \left[\tfrac{1}{2}(\delta_{im}\delta_{jn}+\delta_{in}\delta_{jm}) - \tfrac{1}{3} \delta_{mn}\delta_{ij}\right]
   = 2 s_{mn} - \tfrac{2}{3} s_{kk} \delta_{mn} = 2 s_{mn} \equiv 2 \boldsymbol{s}
$$
</div>
where we have used the observation that $$s_{kk} \equiv \text{tr}(\boldsymbol{s}) = 0$$. We also have
<div>
$$
  \frac{\partial \text{tr}(\boldsymbol{\sigma})}{\partial \boldsymbol{\sigma}} \equiv
  \frac{\partial \sigma_{kk}}{\partial \sigma_{mn}} = \delta_{mn} \equiv \mathbf{I}
$$
</div>

Therefore,
<div>
$$
  \frac{\partial f}{\partial \boldsymbol{\sigma}} =
    \tfrac{1}{\sqrt{2}} \frac{\boldsymbol{s}}{\sqrt{\boldsymbol{s}:\boldsymbol{s}}} -
    \tfrac{1}{3} \frac{dq}{dp}\,\mathbf{I}
  \quad \implies \quad
  \dot{\boldsymbol{\varepsilon}}^p = \dot{\lambda}\,\left[
  \tfrac{1}{\sqrt{2}} \frac{\boldsymbol{s}}{\sqrt{\boldsymbol{s}:\boldsymbol{s}}} -
    \tfrac{1}{3} \frac{dq}{dp}\,\mathbf{I}
  \right]
$$
</div>
{:.notice}

Note that if $$q(p) = 1/\sqrt{3} \sigma_y$$ where $$\sigma_y$$ is the uniaxial stress yield stress, and
there is no hardening, the Drucker-Prager yield condition reduces to the von Mises yield condition.
{:.notice--info}

##### Plane stress Drucker-Prager yield function #####
Recall that in plane stress, the stress tensor can be simplified to the vector
$$\boldsymbol{\sigma} = [\sigma_{11} \, \sigma_{22} \, \sigma_{12}]^T$$.
The mean stress in plane stress can be written as
<div>
$$
  p = \tfrac{1}{3} (\sigma_{11} + \sigma_{22})
$$
</div>
For the deviatoric stress, 
<div>
$$
  \boldsymbol{s} =
    \begin{bmatrix}
      \sigma_{11} & \sigma_{12} & 0 \\ \sigma_{12} & \sigma_{22} & 0 \\ 0 & 0 & 0 
    \end{bmatrix}
    - \tfrac{1}{3} (\sigma_{11}+\sigma_{22})
        \begin{bmatrix}
          1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1
        \end{bmatrix}
$$
</div>
or,
<div>
$$
    \begin{bmatrix}
      s_{11} & s_{12} & 0 \\ s_{12} & s_{22} & 0 \\ 0 & 0 & s_{33} 
    \end{bmatrix} = 
    \tfrac{1}{3}\begin{bmatrix}
       2\sigma_{11} - \sigma_{22} & \sigma_{12} & 0 \\
       \sigma_{12} & -\sigma_{11} + 2\sigma_{22} & 0 \\
       0 & 0 & -\sigma_{11} - \sigma_{22} 
    \end{bmatrix}
$$
</div>
We can write the above in vector notation as
<div>
$$
  \begin{bmatrix}
    s_{11} \\ s_{22} \\ s_{12}
  \end{bmatrix} = \tfrac{1}{3}
  \begin{bmatrix}
    2 & -1 & 0 \\ -1 & 2 & 0 \\ 0 & 0 & 3
  \end{bmatrix}
  \begin{bmatrix}
    \sigma_{11} \\ \sigma_{22} \\ \sigma_{12}
  \end{bmatrix}
  \quad \leftrightarrow \quad
  \boldsymbol{s} = \mathbf{T}\,\boldsymbol{\sigma} \,.
$$
</div>
To express the Drucker-Prager yield function in terms of these vectors, 
using the matrix form of $$\boldsymbol{s}$$, we have
<div>
$$
  J_2 = \tfrac{1}{2} \boldsymbol{s}:\boldsymbol{s} \equiv
      \tfrac{1}{2} \text{tr}(\boldsymbol{s}^T \boldsymbol{s})
    = \tfrac{1}{2} [
          (s_{11}^2 + s_{12}^2) +
          (s_{12}^2 + s_{22}^2) +
          (s_{11} + s_{22})^2 ]
$$
</div>
We can write this in matrix form as
<div>
$$
  \begin{align}
  J_2 & = \tfrac{1}{2}\begin{bmatrix} s_{11} & s_{22} & s_{12} \end{bmatrix}
    \begin{bmatrix} 2 & 1 & 0 \\ 1 & 2 & 0 \\ 0 & 0 & 2 \end{bmatrix}
    \begin{bmatrix} s_{11} \\ s_{22} \\ s_{12} \end{bmatrix} \\
   & = \tfrac{1}{6}\begin{bmatrix} \sigma_{11} & \sigma_{22} & \sigma_{12} \end{bmatrix}
     \begin{bmatrix} 2 & -1 & 0 \\ -1 & 2 & 0 \\ 0 & 0 & 6 \end{bmatrix}
    \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{12} \end{bmatrix}
    =: \tfrac{1}{2} \boldsymbol{\sigma}^T \,\mathbf{P} \,\boldsymbol{\sigma}
  \end{align}
$$
</div>

Therefore, the plane stress Drucker-Prager yield function can be written as
<div>
$$
  f(\boldsymbol{\sigma}) = \sqrt{\tfrac{1}{2} \boldsymbol{\sigma}^T \,\mathbf{P} \,\boldsymbol{\sigma}} - q(p) = 0 
$$
</div>
{:.notice}

This form of the flow yield function is commonly seen in the plasticity literature.

##### Plane stress Drucker-Prager flow rule #####

Recall that the Drucker-Prager flow rule is
<div>
$$
  \dot{\boldsymbol{\varepsilon}}^p = \dot{\lambda}\,\left[
  \tfrac{1}{\sqrt{2}} \frac{\boldsymbol{s}}{\sqrt{\boldsymbol{s}:\boldsymbol{s}}} -
    \tfrac{1}{3} \frac{dq}{dp}\,\mathbf{I}
  \right]
$$
</div>
In matrix form,
<div>
$$
 \begin{bmatrix}
   \dot{\varepsilon}^p_{11} & \dot{\varepsilon}^p_{12} & 0 \\
   \dot{\varepsilon}^p_{12} & \dot{\varepsilon}^p_{22} & 0 \\
   0 & 0 & \dot{\varepsilon}^p_{33} 
 \end{bmatrix} =
 \dot{\lambda} \left(
  \frac{1}{\sqrt{2 \text{tr}(\boldsymbol{s}^T\boldsymbol{s})}}
    \begin{bmatrix}
      s_{11} & s_{12} & 0 \\ s_{12} & s_{22} & 0 \\ 0 & 0 & s_{33}
    \end{bmatrix} - 
    \tfrac{1}{3} \frac{dq}{dp}\,\begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}
   \right)
$$
</div>
We have already seen that
<div>
$$
  \boldsymbol{s}:\boldsymbol{s} =  \text{tr}(\boldsymbol{s}^T\boldsymbol{s})
    = \boldsymbol{\sigma}^T \,\mathbf{P} \,\boldsymbol{\sigma}
$$
</div>
Therefore, in vector form,
<div>
$$
 \begin{bmatrix} \dot{\varepsilon_{11}}^p \\ \dot{\varepsilon_{22}}^p \\ \dot{\varepsilon_{12}}^p \end{bmatrix}
 = \dot{\lambda} \left(
     \frac{1}{\sqrt{2 \boldsymbol{\sigma}^T\,\mathbf{P}\,\boldsymbol{\sigma}}}
     \mathbf{T}\begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{12} \end{bmatrix}
     - \tfrac{1}{3} \frac{dq}{dp}\,\begin{bmatrix} 1 \\ 1 \\ 0 \end{bmatrix}
   \right)
$$
</div>
If we use the standard engineering notation for strain, the flow rule becomes
<div>
$$
 \begin{bmatrix} \dot{\varepsilon_{11}}^p \\ \dot{\varepsilon_{22}}^p \\ 2\dot{\varepsilon_{12}}^p \end{bmatrix}
 = \dot{\lambda} \left(
     \frac{\mathbf{P}\,\boldsymbol{\sigma}}{\sqrt{2 \boldsymbol{\sigma}^T\,\mathbf{P}\,\boldsymbol{\sigma}}}
     - \tfrac{1}{3} \frac{dq}{dp}\,\begin{bmatrix} 1 \\ 1 \\ 0 \end{bmatrix}
   \right)
$$
</div>

Therefore, the plane stress Drucker-Prager flow rule in vector form is
<div>
$$
  \dot{\boldsymbol{\varepsilon}}^p = 
    \dot{\lambda} \left(
     \tfrac{1}{\sqrt{2}}
       \frac{\mathbf{P}\,\boldsymbol{\sigma}}{
         \sqrt{\boldsymbol{\sigma}^T\,\mathbf{P}\,\boldsymbol{\sigma}}
       }
     - \tfrac{1}{3} \frac{dq}{dp}\,\mathbf{I}
   \right)
$$
</div>
{:.notice}

#### Remarks ####
We will explore the return algorithm for the plane stress Drucker-Prager plasticity model in the next
article in this series.

