---
layout: docs
title:  "The plane stress return algorithm"
subheadline: "Biswajit Banerjee"
description: "Part 2 of the Drucker-Prager return algorithm"
date:  2017-03-05 10:30:00
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
In the [previous part]({{ site.url }}/mechanics/plasticity/algorithm/plane-stress-drucker-prager/)
of this discussion, I derived plane stress expressions for linear elasticity, the Drucker-Prager
yield function, and the associated flow rule.

Let us now review the approach used for finding the parameter $$\dot{\lambda}$$ that is needed
in return algorithms, from a purely algebraic standpoint.

Warning: We are dealing only with perfect plasticity (no hardening) in this
series of articles.
{:.notice--warning}

##### Review of 3D plasticity #####
If we assume an additive decomposition of the elastic and plastic strains, the linear elastic
stress-strain relation can be expressed in tensor notation as
<div>
$$
  \boldsymbol{\sigma} = \mathbf{C} : \boldsymbol{\varepsilon}^e
    = \mathbf{C} : (\boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^p)
$$
</div>
The flow rule can be written as
<div>
$$
  \dot{\boldsymbol{\varepsilon}}^p = \dot{\lambda} \frac{\partial f}{\partial\boldsymbol{\sigma}} =: \dot{\lambda} \boldsymbol{n}
$$
</div>
and the consistency condition is
<div>
$$
  \dot{f} = 0 \,.
$$
</div>
From the consistency condition, the flow rule, and the stress-strain relation, we have
<div>
$$
  \dot{f} = \frac{\partial f}{\partial \boldsymbol{\sigma}} : \dot{\boldsymbol{\sigma}}
    = \boldsymbol{n} :
      \mathbf{C} : (\dot{\boldsymbol{\varepsilon}} - \dot{\boldsymbol{\varepsilon}}^p)
    = \boldsymbol{n} :
      \mathbf{C} : (\dot{\boldsymbol{\varepsilon}} - \dot{\lambda} \boldsymbol{n}) = 0
$$
</div>
Solving for $$\dot{\lambda}$$, we have
<div>
$$
  \dot{\lambda} = \frac{\boldsymbol{n} :
      \mathbf{C} : \dot{\boldsymbol{\varepsilon}}}{ \boldsymbol{n} :
      \mathbf{C} : \boldsymbol{n}} \ge 0
$$
</div>
{:.notice}

The elastic-plastic tangent modulus can be computed from the stress-strain relation:
<div>
$$
  \dot{\boldsymbol{\sigma}} = \mathbf{C} : (\dot{\boldsymbol{\varepsilon}} - \dot{\boldsymbol{\varepsilon}}^p) = \mathbf{C} : (\dot{\boldsymbol{\varepsilon}} - \dot{\lambda} \boldsymbol{n})
  =  \mathbf{C} : \left(\dot{\boldsymbol{\varepsilon}} -  \frac{\boldsymbol{n} :
      \mathbf{C} : \dot{\boldsymbol{\varepsilon}}}{ \boldsymbol{n} :
      \mathbf{C} : \boldsymbol{n}} \boldsymbol{n}\right)
$$
</div>
or,
<div>
$$
  \dot{\boldsymbol{\sigma}} = \left[\mathbf{C} -
    \frac{
      (\mathbf{C}: \boldsymbol{n}) \otimes \left(\boldsymbol{n} : \mathbf{C}\right)
    }{
      \boldsymbol{n} : \mathbf{C} : \boldsymbol{n}
    }\right] : \dot{\boldsymbol{\varepsilon}} =: \mathbf{C}^{ep} : \dot{\boldsymbol{\varepsilon}}
$$
</div>
{:.notice}

We will now specialize these results for plane stress conditions.

##### Plane stress plasticity #####
For plane stress Drucker-Prager non-hardening associated plasticity, the above relations
can be expressed in terms of matrices and vectors.

The linear elastic stress-strain relation becomes
<div>
$$
  \boldsymbol{\sigma} = \mathbf{C}\, \boldsymbol{\varepsilon}^e = 
    \mathbf{C}\, (\boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^p)
$$
</div>
or
<div>
$$
  \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{12} \end{bmatrix} =
  \begin{bmatrix} \frac{4\mu(\lambda+\mu)}{\lambda+2\mu} & \frac{2\mu\lambda}{\lambda+2\mu} & 0 \\
                  \frac{2\mu\lambda}{\lambda+2\mu} & \frac{4\mu(\lambda+\mu)}{\lambda+2\mu} & 0 \\
                  0 & 0 & \mu
  \end{bmatrix}
  \begin{bmatrix} \varepsilon_{11} - \varepsilon_{11}^p \\
                  \varepsilon_{22} - \varepsilon_{22}^p \\ 
                  2(\varepsilon_{12} - \varepsilon_{12}^p) \end{bmatrix}
$$
</div>
The yield function is
<div>
$$
  f(\boldsymbol{\sigma}) = \sqrt{\tfrac{1}{2} \boldsymbol{\sigma}^T \,\mathbf{P} \,\boldsymbol{\sigma}} - q(p) = 0 
$$
</div>
and the flow rule is
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
<div>$$\boldsymbol{n} = \begin{bmatrix} n_{11} & n_{22} & n_{12} \end{bmatrix}^T \,.$$</div>
</div>
From the consistency condition, the flow rule, and the stress-strain relation, we have
<div>
$$
  \begin{align}
  \dot{f} & = n_{11}\dot{\sigma}_{11} + n_{12}\dot{\sigma}_{12} +
              n_{21}\dot{\sigma}_{21} + n_{22}\dot{\sigma}_{22}  = 
              \begin{bmatrix} n_{11} & n_{22} & 2n_{12} \end{bmatrix}
              \begin{bmatrix} \dot{\sigma}_{11} \\ \dot{\sigma}_{22} \\ \dot{\sigma}_{12} \end{bmatrix} \\
          & = \begin{bmatrix} n_{11} & n_{22} & 2n_{12} \end{bmatrix}
              \begin{bmatrix} C_{11} & C_{12} & 0 \\ C_{12} & C_{22} & 0 \\ 0 & 0 & C_{33}
              \end{bmatrix}
              \begin{bmatrix} \dot{\varepsilon}_{11} - \dot{\lambda} n_{11} \\
                              \dot{\varepsilon}_{22} - \dot{\lambda} n_{22} \\ 
                              2(\dot{\varepsilon}_{12} - \dot{\lambda} n_{12}) \end{bmatrix}
           = 0
  \end{align}
$$
</div>
Solving for $$\dot{\lambda}$$, we have
<div>
$$
  \dot{\lambda} = \frac{\mathbf{n}^T\,\mathbf{C}\,\dot{\boldsymbol{\varepsilon}}}
                       {\mathbf{n}^T\,\mathbf{C}\,\mathbf{n}} \ge 0 \quad \text{where} \quad
  \mathbf{n}^T = \begin{bmatrix} n_{11} & n_{22} & 2n_{12} \end{bmatrix} \,.
$$
</div>
{:.notice}

The elastic-plastic tangent modulus is computed from the stress-strain relation
<div>
$$
  \dot{\boldsymbol{\sigma}} =
    \mathbf{C}\,(\dot{\boldsymbol{\varepsilon}} - \dot{\boldsymbol{\varepsilon}}^p) =
    \mathbf{C}\,(\dot{\boldsymbol{\varepsilon}} - \dot{\lambda} \boldsymbol{n})
  =  \mathbf{C}\,\left(\dot{\boldsymbol{\varepsilon}} -
       \frac{\mathbf{n}^T\,\mathbf{C}\,\dot{\boldsymbol{\varepsilon}}}
            {\mathbf{n}^T\,\mathbf{C}\,\mathbf{n}} \,\boldsymbol{n}\right)
$$
</div>
to get
<div>
$$
  \dot{\boldsymbol{\sigma}} =
     \left(\mathbf{C} -
          \frac{\mathbf{C}\,\boldsymbol{n}\,\mathbf{n}^T\,\mathbf{C}}
               {\mathbf{n}^T\,\mathbf{C}\,\mathbf{n}}
     \right) \dot{\boldsymbol{\varepsilon}} = \mathbf{C}^{ep} \, \dot{\boldsymbol{\varepsilon}}
$$
</div>
{:.notice}

Vectors $$\boldsymbol{n}$$ and $$\mathbf{n}$$ are related by
{:.notice--warning}
<div>
$$
  \mathbf{n} = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 2 \end{bmatrix} \boldsymbol{n}
$$
</div>
{:.notice--warning}

Next we will see the return algorithm for this problem.

#### Forward Euler Return algorithm ####
Let us define a trial elastic stress rate
<div>
$$
  \dot{\boldsymbol{\sigma}}^{\text{trial}} = \mathbf{C} : \dot{\boldsymbol{\varepsilon}}
$$
</div>
We can express the stress-strain rate relation as
<div>
$$
  \dot{\boldsymbol{\sigma}} = \mathbf{C} : \dot{\boldsymbol{\varepsilon}}^e 
      = \mathbf{C} : (\dot{\boldsymbol{\varepsilon}} - \dot{\boldsymbol{\varepsilon}}^p )
      = \mathbf{C} : (\dot{\boldsymbol{\varepsilon}} - \dot{\lambda}\boldsymbol{n} )
      = \dot{\boldsymbol{\sigma}}^{\text{trial}} - \dot{\lambda}\,\mathbf{C} : \boldsymbol{n}
$$
</div>
From the consistency condition
<div>
$$
  \dot{f} = \boldsymbol{n} : \dot{\boldsymbol{\sigma}}
      = \boldsymbol{n} : (
        \dot{\boldsymbol{\sigma}}^{\text{trial}} - \dot{\lambda}\,\mathbf{C} : \boldsymbol{n})
      = 0
$$
</div>
Therefore,
<div>
$$
  \dot{\lambda} = \frac{\boldsymbol{n} : \dot{\boldsymbol{\sigma}}^{\text{trial}}}
                       {\boldsymbol{n} : \mathbf{C} : \boldsymbol{n}}
$$
</div>
And we get the stress rate in terms of the trial stress rate:
<div>
$$
  \dot{\boldsymbol{\sigma}} = 
      \dot{\boldsymbol{\sigma}}^{\text{trial}} -
      \left(\frac{\boldsymbol{n} : \dot{\boldsymbol{\sigma}}^{\text{trial}}}
                       {\boldsymbol{n} : \mathbf{C} : \boldsymbol{n}}\right)
      \,\mathbf{C} : \boldsymbol{n}
   = \left[\mathbf{I}^{(4)} - \frac{(\boldsymbol{n}:\mathbf{C})\otimes\boldsymbol{n}}
                             {\boldsymbol{n} : \mathbf{C} : \boldsymbol{n}}\right] : 
      \dot{\boldsymbol{\sigma}}^{\text{trial}} 
$$
</div>
{:.notice}
In the discrete form of the return algorithm we assume that we know the state at time $$t_n$$ and
will calculate the state at time $$t_{n+1}$$.  To integrate the above equation we 
go back to matrix notation and observe that
<div>
$$
  \boldsymbol{\sigma}_{n+1}^{\text{trial}} = \boldsymbol{\sigma}_n +
    \mathbf{C}\,(\boldsymbol{\varepsilon}_{n+1} - \boldsymbol{\varepsilon}_{n})
   = \mathbf{C}\,(\boldsymbol{\varepsilon}_{n+1} - \boldsymbol{\varepsilon}_{n}^p)
$$
</div>

The actual stress that we want to compute is
<div>
$$
  \boldsymbol{\sigma}_{n+1} = \mathbf{C}\,\boldsymbol{\varepsilon}_{n+1}^e 
   = \mathbf{C}\,(\boldsymbol{\varepsilon}_{n+1} - \boldsymbol{\varepsilon}_{n+1}^p)
$$
</div>

We integrate the flow rule using forward Euler to get
<div>
$$
  \boldsymbol{\varepsilon}_{n+1}^p =
    \boldsymbol{\varepsilon}_{n}^p + \dot{\lambda} \Delta t \, \boldsymbol{n}_{n} = 
    \boldsymbol{\varepsilon}_{n}^p + \Delta\lambda \, \boldsymbol{n}_{n}
$$
</div>

Therefore, we have
<div>
$$
  \boldsymbol{\sigma}_{n+1} 
   = \mathbf{C}\,(\boldsymbol{\varepsilon}_{n+1} -
    \boldsymbol{\varepsilon}_{n}^p - \Delta\lambda \, \boldsymbol{n}_{n})
   = \boldsymbol{\sigma}_{n+1}^{\text{trial}} - \Delta\lambda \,\mathbf{C}\, \boldsymbol{n}_{n}
$$
</div>
{:.notice}

We don't know $$\Delta\lambda$$ and $$\boldsymbol{\sigma}_{n+1}$$ yet.
To find $$\Delta\lambda$$ we use the consistency condition:
<div>
$$
  f(\boldsymbol{\sigma}_{n+1}) = 0 = 
  \sqrt{\tfrac{1}{2} \boldsymbol{\sigma}_{n+1}^T \,\mathbf{P} \,\boldsymbol{\sigma}_{n+1}} - q(p_{n+1}) 
$$
</div>
At this stage we are ready to get into some unpleasant algebra.

#### Remarks ####
We will explore the rest of the algebra and the backward Euler return algorithm in the next
article in this series after a short break.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

