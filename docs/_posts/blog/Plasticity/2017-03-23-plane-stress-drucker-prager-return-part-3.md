---
title:  "Plane stress return: Spectral decomposition"
subheadline: "Biswajit Banerjee"
description: ""
date:  2017-03-23 10:30:00
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
The previous parts of this series dealt with:

* Part 1: [Background of 3D and plane stress Drucker-Prager plasticity]({{ site.baseurl }}/mechanics/plasticity/algorithm/plane-stress-drucker-prager/)
* Part 2: [A forward Euler return algorithm for plane stress plasticity]({{site.baseurl }}/mechanics/plasticity/algorithm/plane-stress-drucker-prager-return/)
<!--more-->

Let us now take a slight detour and examine the spectral decomposition of the stiffness matrix
($$\mathbf{C}$$) and the transformation matrix ($$\mathbf{P}$$).

##### Spectral decomposition of $$\mathbf{C}$$ #####
Recall that under plane stress conditions the stiffness matrix reduces to
<div>
$$
  \mathbf{C} =
  \begin{bmatrix}
    \frac{4\mu(\lambda+\mu)}{\lambda+2\mu} & \frac{2\mu\lambda}{\lambda+2\mu} & 0 \\
    \frac{2\mu\lambda}{\lambda+2\mu} & \frac{4\mu(\lambda+\mu)}{\lambda+2\mu} & 0 \\
    0 & 0 & \mu
  \end{bmatrix}
$$
</div>
Since $$\mathbf{C}$$ is a real, symmetric, matrix we can find a [spectral decomposition](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix) of
the matrix via a [singular value decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition) of the form
<div>
$$
  \mathbf{C} = \mathbf{Q}_C \mathbf{L}_C \mathbf{Q}_C^T
$$
</div>
where $$\mathbf{Q}_C$$ is an orthogonal matrix whose columns form the basis of the spectral
decomposition and $$\mathbf{L}_C$$ is a diagonal matrix containing the eigenvalues of the
spectral decomposition.

Because we are lazy and prone to algebra errors, let us turn to Mathematica to solve the
singular value decomposition problem.  The script we use is:
<div>
$$
\begin{align}
&\pmb{\text{C11} = 4 \text{mu} (\text{lambda} + \text{mu})/(\text{lambda} + 2 \text{mu})}\\
&\pmb{\text{C22} = \text{C11}}\\
&\pmb{\text{C33} = \text{mu}}\\
&\pmb{\text{C12} = 2 \text{mu}\, \text{lambda}/(\text{lambda} + 2 \text{mu})}\\
&\pmb{\text{CC} = \{\{\text{C11}, \text{C12}, 0\},\{\text{C12}, \text{C22}, 0\},\{0, 0, \text{C33}\}\}}\\
&\pmb{\text{$\$$Assumptions} =\text{  }\text{lambda}\in  \text{Reals} \,\&\&\,\text{  }\text{mu} \in  \text{Reals} \,\&\& \,\text{mu} > 0 \,\&\&\,}\\
& \quad \pmb{\text{CC} \in \text{Matrices}[\{3,3\},\text{Reals},\text{Symmetric}]}\\
&\pmb{\{\text{QC}, \text{LC}, \text{QCT}\} = \text{SingularValueDecomposition}[\text{CC}]\text{//} \,\text{FullSimplify}} \\
&\pmb{\text{QC} \text{/.} \{\text{lambda} \to  \text{Ee}\, \text{nu}/((1+ \text{nu}) (1 - 2 \text{nu})), \text{mu} \to  \text{Ee}/(2 (1 + \text{nu}))\}\text{//} \text{FullSimplify}} \\
& \pmb{\text{LC} \text{/.} \{\text{lambda} \to  \text{Ee}\, \text{nu}/((1+ \text{nu}) (1 - 2 \text{nu})), \text{mu} \to  \text{Ee}/(2 (1 + \text{nu}))\}\text{//} \text{FullSimplify}}
\end{align}
$$
</div>
The result produced by Mathematica is
<div>
$$
\begin{align}
 \mathbf{Q}_C &= \left\{\left\{0,\frac{\text{Sign}[1+\text{nu}]}{\sqrt{2} \text{Sign}[1-\text{nu}]},-\frac{1}{\sqrt{2}}\right\},\left\{0,\frac{\text{Sign}[1+\text{nu}]}{\sqrt{2}
\text{Sign}[1-\text{nu}]},\frac{1}{\sqrt{2}}\right\},\{1,0,0\}\right\} \\
 \mathbf{L}_C & = \left\{\left\{\frac{\text{Ee}}{2+2 \text{nu}},0,0\right\},\left\{0,\frac{\text{Ee} \text{Abs}\left[\frac{1+\text{nu}}{1-\text{nu}}\right]}{1+\text{nu}},0\right\},\left\{0,0,\frac{\text{Ee}}{1+\text{nu}}\right\}\right\}
\end{align}
$$
</div>
Since the Poisson's ratio $$\nu \in [-1, 1/2]$$, we have $$\text{sign}(1 + \nu) = 1$$ and $$\text{sign}(1 - \nu) = 1$$.  Therefore,
<div>
 $$
   \mathbf{Q}_C = \begin{bmatrix} 0 & \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}} \\
                         0 & \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} \\
                         1 & 0 & 0 \end{bmatrix}
   \quad
   \mathbf{L}_C = \begin{bmatrix} \frac{E}{2(1+\nu)} & 0 & 0 \\ 0 & \frac{E}{1-\nu} & 0 \\
                                  0 & 0 & \frac{2E}{2(1+\nu)} \end{bmatrix}
    = \begin{bmatrix} \mu & 0 & 0 \\ 0 & \frac{E}{1-\nu} & 0 \\
                                  0 & 0 & 2\mu \end{bmatrix}
 $$
</div>

##### Spectral decomposition of $$\mathbf{P}$$ #####
The transformation matrix $$\mathbf{P}$$ has components
<div>
$$
  \mathbf{P} = \begin{bmatrix} 2 & -1 & 0 \\ -1 & 2 & 0 \\ 0 & 0 & 6 \end{bmatrix}
$$
</div>
As we did for $$\mathbf{C}$$, we can perform a spectral decomposition of
matrix $$\mathbf{P}$$ into
<div>
$$
  \mathbf{P} = \mathbf{Q}_P \mathbf{L}_P \mathbf{Q}_P^T
$$
</div>
Proceeding as before with Mathematica, we find that
<div>
$$
   \mathbf{Q}_P = \begin{bmatrix} 0 & \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}} \\
                         0 & \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} \\
                         1 & 0 & 0 \end{bmatrix}
   \quad
   \mathbf{L}_P = \begin{bmatrix} 6 & 0 & 0 \\ 0 & 1 & 0 \\
                                  0 & 0 & 3 \end{bmatrix}
$$
</div>

##### Product of $$\mathbf{P}$$ and $$\mathbf{C}$$ #####
From the above we see that $$\mathbf{Q}_C = \mathbf{Q}_P =: \mathbf{Q}$$.  Therefore,
<div>
$$
  \mathbf{P}\mathbf{C} = \mathbf{Q} \mathbf{L}_P \mathbf{Q}^T  \mathbf{Q} \mathbf{L}_C \mathbf{Q}^T
     = \mathbf{Q} \mathbf{L}_P \mathbf{L}_C \mathbf{Q}^T
     = \mathbf{Q} \mathbf{L}_C \mathbf{L}_P \mathbf{Q}^T
     = \mathbf{C} \mathbf{P}
$$
</div>

#### Remarks ####
We will use these decompositions to diagonalize the complicated matrix expressions in the
next part of this series.


