---
title:  "Plane stress forward Euler Drucker-Prager"
subheadline: "Biswajit Banerjee"
description: ""
date:  2017-03-24 10:30:00
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
In [Part 2]({{site.baseurl }}/mechanics/plasticity/algorithm/plane-stress-drucker-prager-return/)
of this series, we saw that a forward Euler return algorithm leads to the discretized equations
<!--more-->
<div>
$$
 \boldsymbol{\sigma}_{n+1}
   = \boldsymbol{\sigma}_{n+1}^{\text{trial}} - \Delta\lambda \,\mathbf{C}\, \boldsymbol{n}_{n}
$$
</div>
{:.notice}
where
<div>
$$
 \boldsymbol{n}_n = \begin{bmatrix} n^n_{11} \\ n^n_{22} \\ n^n_{12} \end{bmatrix}
   = \left(
     \tfrac{1}{\sqrt{2}}
       \frac{\mathbf{P}\,\boldsymbol{\sigma_n}}{
         \sqrt{\boldsymbol{\sigma_n}^T\,\mathbf{P}\,\boldsymbol{\sigma_n}}
       }
     - \tfrac{1}{3} \frac{dq}{dp}(\boldsymbol{\sigma}_n)\,\mathbf{I}
   \right)
$$
</div>
and
<div>
$$
  \mathbf{C} = \begin{bmatrix}
    \frac{4\mu(\lambda+\mu)}{\lambda+2\mu} & \frac{2\mu\lambda}{\lambda+2\mu} & 0 \\
    \frac{2\mu\lambda}{\lambda+2\mu} & \frac{4\mu(\lambda+\mu)}{\lambda+2\mu} & 0 \\
    0 & 0 & \mu
  \end{bmatrix}
  \quad \mathbf{P} = \begin{bmatrix}
                       2 & -1 & 0 \\ -1 & 2 & 0 \\ 0 & 0 & 6
                      \end{bmatrix}
   \quad \boldsymbol{\sigma}_n = 
     \begin{bmatrix} \sigma^n_{11} \\ \sigma^n_{22} \\ \sigma^n_{12} \end{bmatrix}
   \quad \mathbf{I} =  \begin{bmatrix} 1 \\ 1 \\ 0 \end{bmatrix}
$$
</div>
The discrete consistency condition is
<div>
$$
  f(\boldsymbol{\sigma}_{n+1}) = 0 = 
  \sqrt{\tfrac{1}{2} \boldsymbol{\sigma}_{n+1}^T \,\mathbf{P} \,\boldsymbol{\sigma}_{n+1}} - q(p_{n+1}) 
$$
</div>
{:.notice}
In [Part 3]({{ site.baseurl }}/mechanics/plasticity/algorithm/plane-stress-drucker-prager-return-part-3/), we saw that
<div>
$$
  \mathbf{C} = \mathbf{Q} \mathbf{L}_C \mathbf{Q}^T \quad
  \mathbf{P} = \mathbf{Q} \mathbf{L}_P \mathbf{Q}^T 
$$
</div>
{:.notice}
where
<div>
$$
  \mathbf{Q} = \begin{bmatrix} 0 & \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}} \\
                         0 & \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} \\
                         1 & 0 & 0 \end{bmatrix}
  \quad
  \mathbf{L}_P = \begin{bmatrix} 6 & 0 & 0 \\ 0 & 1 & 0 \\
                                  0 & 0 & 3 \end{bmatrix}
  \quad
  \mathbf{L}_C = \begin{bmatrix} \mu & 0 & 0 \\ 0 & \frac{E}{1-\nu} & 0 \\
                                  0 & 0 & 2\mu \end{bmatrix}
$$
</div>
Let us now try to find $$\Delta\lambda$$ and $$\boldsymbol{\sigma}_{n+1}$$ and check
whether we can find simple expressions for those.

##### Finding $$\Delta\lambda$$ #####
To find $$\Delta\lambda$$, we use the consistency condition and express all quantities in
that expression in terms of the trial stress.  Then,
<div>
$$
  \tfrac{1}{\sqrt{2}}
  \sqrt{
  \left[\boldsymbol{\sigma}^\text{trial}_{n+1} -\Delta\lambda \mathbf{C}\boldsymbol{n}_n\right]^T
  \mathbf{P}
  \left[\boldsymbol{\sigma}^\text{trial}_{n+1} -\Delta\lambda \mathbf{C}\boldsymbol{n}_n\right]}
  = q(\text{tr}[\boldsymbol{\sigma}_{n+1}^\text{trial} - \Delta\lambda\mathbf{C}\boldsymbol{n}_n])
$$
</div>
or
<div>
$$
  \left[(\boldsymbol{\sigma}^\text{trial}_{n+1})^T -\Delta\lambda \boldsymbol{n}_n^T\mathbf{C}\right]
  \left[\mathbf{P}\boldsymbol{\sigma}^\text{trial}_{n+1} -\Delta\lambda \mathbf{P}\mathbf{C}\boldsymbol{n}_n\right]
  = 2q^2(\text{tr}[\boldsymbol{\sigma}_{n+1}^\text{trial} - \Delta\lambda\mathbf{C}\boldsymbol{n}_n])
$$
</div>
or
<div>
$$
  (\boldsymbol{\sigma}^\text{trial}_{n+1})^T\mathbf{P}\boldsymbol{\sigma}^\text{trial}_{n+1}
  - 2\Delta\lambda (\boldsymbol{\sigma}^\text{trial}_{n+1})^T\mathbf{P}\mathbf{C}\boldsymbol{n}_n
  + (\Delta\lambda)^2 (\mathbf{C}\boldsymbol{n}_n)^T\mathbf{P}\mathbf{C}\boldsymbol{n}_n \\
   = 2q^2(\text{tr}[\boldsymbol{\sigma}_{n+1}^\text{trial} - \Delta\lambda\mathbf{C}\boldsymbol{n}_n])
$$
</div>
{:.notice}

Because of the potential dependence on $$\Delta\lambda$$ in the expression for $$q^2$$,
this equation is best solved
using a root-finding method.  The quadratic equation in $$\Delta\lambda$$ has two solutions, only
one of which should be greater than 0 for the solution to be unique.  It is not straightforward
to show that this is indeed the case.
{:.notice--warning}


##### An attempt at simplification #####
Let us try to find a simpler expression for $$\mathbf{C}\boldsymbol{n}_n$$ and related
quantities.  We have
<div>
$$
  \mathbf{C}\boldsymbol{n}_n = 
    \left(
     \tfrac{1}{\sqrt{2}}
       \frac{\mathbf{C}\,\mathbf{P}\,\boldsymbol{\sigma_n}}{
         \sqrt{\boldsymbol{\sigma_n}^T\,\mathbf{P}\,\boldsymbol{\sigma_n}}
       }
     - \tfrac{1}{3} \frac{dq}{dp}(\boldsymbol{\sigma}_n)\,\mathbf{C}\,\mathbf{I}
   \right)
$$
</div>
Therefore,
<div>
$$
  \mathbf{P}\mathbf{C}\boldsymbol{n}_n = 
    \left(
     \tfrac{1}{\sqrt{2}}
       \frac{\mathbf{P}\mathbf{C}\,\mathbf{P}\,\boldsymbol{\sigma_n}}{
         \sqrt{\boldsymbol{\sigma_n}^T\,\mathbf{P}\,\boldsymbol{\sigma_n}}
       }
     - \tfrac{1}{3} \frac{dq}{dp}(\boldsymbol{\sigma}_n)\,\mathbf{P}\mathbf{C}\,\mathbf{I}
   \right)
$$
</div>
and
<div>
$$
  (\mathbf{C}\boldsymbol{n}_n)^T\mathbf{P}\mathbf{C}\boldsymbol{n}_n = 
    \left(
     \tfrac{1}{\sqrt{2}}
       \frac{\mathbf{C}\,\mathbf{P}\,\boldsymbol{\sigma_n}}{
         \sqrt{\boldsymbol{\sigma_n}^T\,\mathbf{P}\,\boldsymbol{\sigma_n}}
       }
     - \tfrac{1}{3} \frac{dq}{dp}(\boldsymbol{\sigma}_n)\,\mathbf{C}\,\mathbf{I}
   \right)^T
    \left(
     \tfrac{1}{\sqrt{2}}
       \frac{\mathbf{P}\mathbf{C}\,\mathbf{P}\,\boldsymbol{\sigma_n}}{
         \sqrt{\boldsymbol{\sigma_n}^T\,\mathbf{P}\,\boldsymbol{\sigma_n}}
       }
     - \tfrac{1}{3} \frac{dq}{dp}(\boldsymbol{\sigma}_n)\,\mathbf{P}\mathbf{C}\,\mathbf{I}
   \right)
$$
</div>
Now,
<div>
$$
  \mathbf{C}\mathbf{P}\boldsymbol{\sigma}_n = \mathbf{Q}\mathbf{L}_C\mathbf{L}_P\mathbf{Q}^T\boldsymbol{\sigma}_n \quad \text{and} \quad
  \boldsymbol{\sigma}_n^T\mathbf{P}\boldsymbol{\sigma}_n = \boldsymbol{\sigma}_n^T\mathbf{Q}\mathbf{L}_P\mathbf{Q}^T\boldsymbol{\sigma}_n 
$$
</div>
Define $$\mathbf{R}_n := \mathbf{Q}^T\boldsymbol{\sigma}_n$$ and
$$\mathbf{R}_{n+1}^{\text{trial}} := \mathbf{Q}^T\boldsymbol{\sigma}_{n+1}^{\text{trial}}$$.
Then,
<div>
$$
  \mathbf{C}\mathbf{P}\boldsymbol{\sigma}_n = \mathbf{Q}\mathbf{L}_C\mathbf{L}_P\mathbf{R}_n
  \quad , \quad
  \mathbf{P}\mathbf{C}\mathbf{P}\boldsymbol{\sigma}_n = \mathbf{Q}\mathbf{L}_C\mathbf{L}_P^2\mathbf{R}_n
  \quad \text{and} \quad
  \boldsymbol{\sigma}_n^T\mathbf{P}\boldsymbol{\sigma}_n = \mathbf{R}_n^T\mathbf{L}_P\mathbf{R}_n
$$
</div>
Then we get the following expressions for the terms in the nonlinear equation for
$$\Delta\lambda$$ (after defining $$\mathbf{Q}_I := \mathbf{Q}^T \mathbf{I}$$):
<div>
$$
 \begin{align}
  & (\boldsymbol{\sigma}^\text{trial}_{n+1})^T\mathbf{P}\boldsymbol{\sigma}^\text{trial}_{n+1}
  = (\mathbf{R}_{n+1}^\text{trial})^T\mathbf{L}_P\mathbf{R}_{n+1}^\text{trial} \\
  & (\boldsymbol{\sigma}^\text{trial}_{n+1})^T\mathbf{P}\mathbf{C}\boldsymbol{n}_n
  = (\mathbf{R}_{n+1}^{\text{trial}})^T\mathbf{L}_C\mathbf{L}_P
    \left(
     \tfrac{1}{\sqrt{2}}
       \frac{\mathbf{L}_P\mathbf{R}_n}{\sqrt{\mathbf{R}_n^T\mathbf{L}_P\mathbf{R}_n}}
     - \tfrac{1}{3} \frac{dq}{dp}\,\mathbf{Q}_I
   \right)\\
  & (\mathbf{C}\boldsymbol{n}_n)^T\mathbf{P}\mathbf{C}\boldsymbol{n}_n = 
  \tfrac{1}{2} \frac{\mathbf{R}_n^T\mathbf{L}_C^2\mathbf{L}_P^3\mathbf{R}_n}{\sqrt{\mathbf{R}_n^T\mathbf{L}_P\mathbf{R}_n}} -
  \tfrac{2}{3\sqrt{2}} \frac{dq}{dp}\frac{\mathbf{R}_n^T\mathbf{L}_C^2\mathbf{L}_P^2}{\sqrt{\mathbf{R}_n^T\mathbf{L}_P\mathbf{R}_n}}\mathbf{Q}_{I} +
  \tfrac{1}{9}\left(\frac{dq}{dp}\right)^2 \mathbf{Q}_I^T\mathbf{L}_C^2\mathbf{L}_P\mathbf{Q}_I 
 \end{align}
$$
</div>
{:.notice}

For $$J_2$$ (von Mises) plasticity, the $$dq/dp$$ terms are zero and further simplification
is possible.  However, for a general Drucker-Prager material no such benefit is found from
the use of the spectral decomposition, other than efficiency of computation.

In our case, the two column matrices that can help with computational efficiency are
$$\mathbf{R}$$ and $$\mathbf{Q}_I$$:
<div>
$$
  \mathbf{R} = \tfrac{1}{\sqrt{2}}\begin{bmatrix} \sqrt{2}\sigma_{12} \\ \sigma_{11}+\sigma_{22}
    \\ -\sigma_{11}+\sigma_{22} \end{bmatrix}
  \quad \text{and} \quad
  \mathbf{Q}_I = \sqrt{2} \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix} \,.
$$
</div>
{:.notice--info}

#### Remarks ####
These results show that the "radial return" approach is applicable only to a very
specific class of plasticity models.  The next parts of this series will discuss
the backward Euler approach and geometric approaches for return algorithms.

