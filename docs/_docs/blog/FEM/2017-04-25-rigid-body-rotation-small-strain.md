---
layout: docs
title:  "Are stresses tensile or compressive during rigid body rotation?"
subheadline: "Biswajit Banerjee"
description: "Small strain finite elements"
date:  2017-04-25 10:30:00
categories:
    - FEM
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

#### Introduction ####
Recently, I came across a question in StackExchange that pointed out that some
books on continuum mechanics suggest that an element will *increase in size*
when rotated if a small strain approximation is used in a finite element simulation.
This issue may be a source of confusion for students of mechanics.  Let us
explore some aspects of the problem and try to answer the StackExchange question.  

#### The question ####
The animation below shows the motion of a two-dimensional square in the $$x$$-$$y$$-plane.
We rotate the square clockwise around one of its corners and the angle of rotation is $$\beta$$.
If the position of a point in the reference configuration of the square is $$(X,Y)$$,
the rotated position $$(x,y)$$ is
<div>
$$
  \begin{bmatrix} x \\ y \end{bmatrix} = \begin{bmatrix} \cos\beta & \sin\beta \\ -\sin\beta & \cos\beta \end{bmatrix} \begin{bmatrix} X \\ Y \end{bmatrix}
$$
</div>
Using the definition of engineering strains: 
<div>
$$
  \varepsilon_{ij} = \tfrac{1}{2}\left(\frac{\partial u_i}{\partial X_j} + \frac{\partial u_j}{\partial X_i}\right)
$$
</div>
and the displacements
<div>
$$
  \begin{bmatrix} u_x \\ u_y \end{bmatrix} = \begin{bmatrix} x - X \\ y - Y \end{bmatrix}
    = \begin{bmatrix} X\cos\beta + Y\sin\beta - X \\ -X\sin\beta + Y\cos\beta - Y\end{bmatrix}
$$
</div>
we have
<div>
$$
  \begin{bmatrix} \varepsilon_{xx} \\ \varepsilon_{yy} \\ 2\varepsilon_{xy} \end{bmatrix} =
  \begin{bmatrix} \cos\beta - 1 \\ \cos\beta - 1 \\ 0 \end{bmatrix} 
$$
</div>
The animation below shows the implications of the above.  The blue square is the
reference configuration.  Pure rotation of the square leads to the configuration in green.
The pink square shows the implication of the strains calculated above on the shape of
the square.
<div>
  <input name="restartButton" type="button" value="Restart animation" onclick="restartSimulation()" />
</div>
<div>
  <canvas id="rot-elem" height="450" width="450"></canvas>
</div>

<div>
To depict the deformed shape with the apparent strains we have used the transformations
$$
  \begin{bmatrix} \varepsilon_{xx}^{\text{rot}} \\ \varepsilon_{yy}^{\text{rot}} \\ 2\varepsilon_{xy}^{\text{rot}} \end{bmatrix} =
  \begin{bmatrix} \cos^2\beta & \sin^2\beta & \sin\beta\cos\beta \\
                  \sin^2\beta & \cos^2\beta & -\sin\beta\cos\beta \\
                  -2\sin\beta\cos\beta & 2\sin\beta\cos\beta & \cos^2\beta-\sin^2\beta
  \end{bmatrix} 
  \begin{bmatrix} \varepsilon_{xx} \\ \varepsilon_{yy} \\ 2\varepsilon_{xy} \end{bmatrix} 
$$
to compute the strains in a coordinate system aligned with the sides of the rotated square.
You can easily check that these strains are identical to those in the unrotated coordinate
system.
</div>
{:.notice--info}

The animation clearly shows, if we consider the strained square, that the square appears
to **decrease** in size (compressive stresses), invert, and then increase in size as the
rotation proceeds.  Some textbooks depict an **increase** (tensile stresses) in the size of
the square to illustrate the same issue, apparently the result of finite element
computations.  The StackExchange question wanted to find out
which was correct, an increase or a decrease.
{:.notice--warning}

#### The finite element solution ####
The question seems to be "why does FEA seemingly produce results that are contradictory 
to theory?"  Let us do the FE calculation for a rotating square and see what the results are.

Consider a square element with nodes 

~~~ bash
id  x    y
1 -1.0 -1.0
2  1.0 -1.0
3  1.0  1.0
4 -1.0  1.0
~~~

The displacement field in the element is
<div>
$$
  u_x(x,y) = \sum_{j=1}^4 u_x^j N_j(x,y) ~,~~
  u_y(x,y) = \sum_{j=1}^4 u_y^j N_j(x,y)
$$
</div>
and the corresponding (small) strain field is
<div>
$$
  \begin{align}
  \varepsilon_{xx}(x,y) & = \sum_{j=1}^4 u_x^j \frac{\partial N_j(x,y)}{\partial x} ~,~~
  \varepsilon_{yy}(x,y) = \sum_{j=1}^4 u_y^j \frac{\partial N_j(x,y)}{\partial y} \\
  \varepsilon_{xy}(x,y) & = \tfrac{1}{2}\sum_{j=1}^4 \left[u_x^j \frac{\partial N_j(x,y)}{\partial y} + u_y^j \frac{\partial N_j(x,y)}{\partial x}\right]
  \end{align}
$$
</div>
The nodal shape functions are
<div>
$$
  \begin{align}
  N_1(x,y) & = \frac{(1-x)(1-y)}{4} ~,~~
  N_2(x,y) = \frac{(1+x)(1-y)}{4} \\
  N_3(x,y) & = \frac{(1+x)(1+y)}{4} ~,~~
  N_4(x,y) = \frac{(1-x)(1+y)}{4}
  \end{align}
$$
</div>
The gradients of the shape functions are
<div>
$$
  \begin{align}
   G_{1x} := \frac{\partial N_1(x,y)}{\partial x} & = -\frac{(1-y)}{4} ~,~~
   G_{1y} := \frac{\partial N_1(x,y)}{\partial y} = -\frac{(1-x)}{4} \\
   G_{2x} := \frac{\partial N_2(x,y)}{\partial x} & = \frac{(1-y)}{4} ~,~~
   G_{2y} := \frac{\partial N_2(x,y)}{\partial y} = -\frac{(1+x)}{4} \\
   G_{3x} := \frac{\partial N_3(x,y)}{\partial x} & = \frac{(1+y)}{4} ~,~~
   G_{3y} := \frac{\partial N_3(x,y)}{\partial y} = \frac{(1+x)}{4} \\
   G_{4x} := \frac{\partial N_4(x,y)}{\partial x} & = -\frac{(1+y)}{4} ~,~~
   G_{4y} := \frac{\partial N_4(x,y)}{\partial y} = \frac{(1-x)}{4}
  \end{align}
$$
</div>
Therefore the strain field can be expressed as
<div>
$$
  \begin{bmatrix}  
   \varepsilon_{xx}(x,y) \\ \varepsilon_{yy}(x,y) \\ 2\varepsilon_{xy}(x,y) 
  \end{bmatrix}
 = \begin{bmatrix} 
     G_{1x} & G_{2x} & G_{3x} & G_{4x} & 0 & 0 & 0 & 0 \\
     0 & 0 & 0 & 0 & G_{1y} & G_{2y} & G_{3y} & G_{4y} \\
     G_{1y} & G_{2y} & G_{3y} & G_{4y} & G_{1x} & G_{2x} & G_{3x} & G_{4x}
   \end{bmatrix} \mathbf{u}
 $$
</div>
 where
<div>
 $$
   \mathbf{u} = \begin{bmatrix}
     u_x^1 & u_x^2 & u_x^3 & u_x^4 & u_y^1 & u_y^2 & u_y^3 & u_y^4
   \end{bmatrix}^T
 $$
</div>
Now consider the situation where the square is rotated by 90 degrees clockwise around
node 1 so that the new positions of the nodes are

~~~ bash
id  x    y
1 -1.0 -1.0
2 -1.0 -3.0
3  1.0 -3.0
4  1.0 -1.0
~~~

*Note that we are not assuming any deformation of the square.*
Then the nodal displacements are

~~~ bash
id  u_x  u_y
1   0.0  0.0
2  -2.0 -2.0
3   0.0 -4.0
4   2.0 -2.0
~~~

Let us do the calculation using Matlab/Octave.
{% highlight matlab %}
>> x = [-1 1 1 -1];
>> y = [-1 -1 1 1];
>> G1x = -(1-y)/4
G1x = -0.50000  -0.50000  -0.00000  -0.00000
>> G1y = -(1-x)/4
G1y = -0.50000  -0.00000  -0.00000  -0.50000
>> G2x = (1-y)/4
G2x = 0.50000   0.50000   0.00000   0.00000
>> G2y = -(1+x)/4
G2y = -0.00000  -0.50000  -0.50000  -0.00000
>> G3x = (1+y)/4
G3x = 0.00000   0.00000   0.50000   0.50000
>> G3y = (1+x)/4
G3y = 0.00000   0.50000   0.50000   0.00000
>> G4x = -(1+y)/4
G4x = -0.00000  -0.00000  -0.50000  -0.50000
>> G4y = (1-x)/4
G4y = 0.50000   0.00000   0.00000   0.50000
>> u = [0 -2 0 2 0 -2 -4 -2];
>> B1 = [[G1x(1) G2x(1) G3x(1) G4x(1) 0 0 0 0];...
         [0 0 0 0 G1y(1) G2y(1) G3y(1) G4y(1)];...
         [G1y(1) G2y(1) G3y(1) G4y(1) G1x(1) G2x(1) G3x(1) G4x(1)]]
B1 =
  -0.50000   0.50000   0.00000  -0.00000   0.00000   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000   0.00000  -0.50000  -0.00000   0.00000   0.50000
  -0.50000  -0.00000   0.00000   0.50000  -0.50000   0.50000   0.00000  -0.00000
>> B2 = [[G1x(2) G2x(2) G3x(2) G4x(2) 0 0 0 0];...
         [0 0 0 0 G1y(2) G2y(2) G3y(2) G4y(2)];...
         [G1y(2) G2y(2) G3y(2) G4y(2) G1x(2) G2x(2) G3x(2) G4x(2)]]
B2 =
  -0.50000   0.50000   0.00000  -0.00000   0.00000   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000   0.00000  -0.00000  -0.50000   0.50000   0.00000
  -0.00000  -0.50000   0.50000   0.00000  -0.50000   0.50000   0.00000  -0.00000
>> B3 = [[G1x(3) G2x(3) G3x(3) G4x(3) 0 0 0 0];...
         [0 0 0 0 G1y(3) G2y(3) G3y(3) G4y(3)];...
         [G1y(3) G2y(3) G3y(3) G4y(3) G1x(3) G2x(3) G3x(3) G4x(3)]]
B3 =
  -0.00000   0.00000   0.50000  -0.50000   0.00000   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000   0.00000  -0.00000  -0.50000   0.50000   0.00000
  -0.00000  -0.50000   0.50000   0.00000  -0.00000   0.00000   0.50000  -0.50000
>> B4 = [[G1x(4) G2x(4) G3x(4) G4x(4) 0 0 0 0];...
         [0 0 0 0 G1y(4) G2y(4) G3y(4) G4y(4)];...
         [G1y(4) G2y(4) G3y(4) G4y(4) G1x(4) G2x(4) G3x(4) G4x(4)]]
B4 =
  -0.00000   0.00000   0.50000  -0.50000   0.00000   0.00000   0.00000   0.00000
   0.00000   0.00000   0.00000   0.00000  -0.50000  -0.00000   0.00000   0.50000
  -0.50000  -0.00000   0.00000   0.50000  -0.00000   0.00000   0.50000  -0.50000
>> eps1 = B1*u'
eps1 =
  -1
  -1
   0
>> eps2 = B2*u'
eps2 =
  -1
  -1
   0
>> eps3 = B3*u'
eps3 =
  -1
  -1
   0
>> eps4 = B4*u'
eps4 =
  -1
  -1
   0
{% endhighlight %}
In the above script, we compute the strains at the four nodes of the square
using the nodal displacement vector

~~~ bash
u = [0 -2 0 2 0 -2 -4 -2]
~~~

and find that the strain at each of the nodes is  

~~~ bash
eps = B*u' = [-1 -1 0]'
~~~~
    
You can easily show that **the finite element solution is identical to the analytical
solution**.   The solution implies that **stresses will develop in the element due to
pure rigid body rotation even if the element does not deform** if we use small strain theory.
{:.notice--info}

This result also implies that the element cannot grow in size with rigid body rotation,
and that the illustrations shown in some textbooks are wrong even though the stress may
transition from more compressive to less compressive during the rotation process.
{:.notice--warning}

#### Remarks ####
The finite element method is typically implemented in displacement form. So one specifies
the displacements and then find the stresses in the element.  An alternative would be to
specify moments that would rotate the element based on forces computed using a displacement
based solution.  

Of course one should not miss the point of the exercise which is to show that the
correct strain measures have to be used for large deformation problems.  Even if a small
strain measure is used, rotations have to be removed from the deformation before any
stress computations are done.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

<script src="{{ site.url }}/assets/js/d3.v4.min.js"></script>
<script src="{{ site.url }}/assets/js/elemRotate.js"></script>
