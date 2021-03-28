---
title:  "Geometric closest point return algorithm"
subheadline: "Biswajit Banerjee"
description: "Part 8 of the series on plasticity return algorithms"
date:  2017-03-31 10:30:00
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
In [Part 7]({{site.baseurl }}/mechanics/plasticity/algorithm/closest-point-return/), we saw that
for isotropic elastic materials and perfect associated plasticity, the trial stress
and the actual stress are at the shortest distance from each other in a transformed stress
space.  We also saw that the transformed stress can be expressed as
<!--more-->
<div>
$$
  \boldsymbol{\sigma}^\star 
  = \frac{z}{\sqrt{3\kappa}}\,\mathbf{E}_z +
    \frac{r}{\sqrt{2\mu}}\,\mathbf{E}_r
$$
</div>
where
<div>
$$
  z = \tfrac{1}{\sqrt{3}}\,\text{tr}(\boldsymbol{\sigma}) ~,~~
  r = \lVert \boldsymbol{s} \rVert 
$$
</div>
and
<div>
$$
  \mathbf{E}_z = \tfrac{1}{\sqrt{3}}\,\mathbf{I} ~,~~
  \mathbf{E}_r = \frac{\boldsymbol{s}}{\lVert\boldsymbol{s}\rVert} \,.
$$
</div>
We can show that the transformed stress vector remains geometrically unchanged (in the sense that
angles are unchanged) if we express it as
<div>
$$
  \boldsymbol{\sigma}^\star 
  = z\,\mathbf{E}_z + \sqrt{\frac{3\kappa}{2\mu}}\,r\,\mathbf{E}_r
  =: z\,\mathbf{E}_z + r'\,\mathbf{E}_r
$$
</div>
{:.notice}
This observation can be used to perform a purely geometric stress update that can
be quite efficient for nonlinear Drucker-Prager plasticity and many other phenomenological
plasticity models.  Let us see how this concept has been applied to the [Arena plasticity
model](https://github.com/bbanerjee/ParSim/blob/master/Vaango/src/CCA/Components/MPM/ConstitutiveModel/Arenisca3PartiallySaturated.cc).  

##### The Arena yield function #####
The Arena model is an extension to partially saturated soils of the Arenisca model
developed by R. M. Brannon and co-workers.  The yield function used by this model is a
nonlinear Drucker-Prager model with a compression cap.  Further details can be found in
the [Arena theory manual](https://github.com/bbanerjee/ParSim/tree/master/Vaango/Manuals/TheoryGuide/ArenaSoil).

If the volumetric and deviatoric components of the total stress are
<div>
$$
    \bar{p} := -\tfrac{1}{3} \text{tr}(\boldsymbol{\sigma}) \quad \text{and} \quad \boldsymbol{s} := \boldsymbol{\sigma} + \bar{p} \mathbf{I}
$$
</div>
we can define
<div>
$$
  \begin{align}
    \bar{p}_\text{eff} & := -\tfrac{1}{3} \text{tr}(\boldsymbol{\sigma}_\text{eff}) = \bar{p} - B \bar{p}^w \\
    \boldsymbol{s}_\text{eff} & := \boldsymbol{\sigma}_\text{eff} + \bar{p}_\text{eff} \mathbf{I} = \boldsymbol{\sigma} + \bar{p} \mathbf{I} = \boldsymbol{s} \\
    J_2^\text{eff} & := \tfrac{1}{2} \boldsymbol{s}_\text{eff}:\boldsymbol{s}_\text{eff}  \,.
  \end{align}
$$
</div>
Then the Arena yield function is
<div>
$$
     f(\boldsymbol{\sigma}, B, \bar{p}^w, \bar{X}) =
       \sqrt{J_2^\text{eff}} - F_f(\bar{p}_\text{eff}) \, F_c(\bar{p}_\text{eff}, \bar{X})
$$
</div>
{:.notice}
where
<div>
$$
    F_f(\bar{p}_\text{eff})  = a_1 - a_3 \exp[- 3 a_2 \bar{p}_\text{eff}] + 3 a_4 \bar{p}_\text{eff}
$$
</div>
{:.notice}
and
<div>
$$
    F_c(\bar{p}_\text{eff}, \bar{X})  =
       \begin{cases}
         1 & \quad \text{for}\quad 3\bar{p}_\text{eff} \le \bar{\kappa} \\
         \sqrt{1 - \left(\cfrac{3\bar{p}_\text{eff} - \bar{\kappa}}{\bar{X}_\text{eff} - \bar{\kappa}}\right)^2} &
           \quad \text{for}\quad 3\bar{p}_\text{eff} > \bar{\kappa} \,.
       \end{cases}
$$
</div>
{:.notice}
Here $$a_i$$ are material parameters, $$\bar{X}_\text{eff}(\boldsymbol{\varepsilon}^p, B, \bar{p}^w) = \bar{X} - 3B\bar{p}^w $$ is the
shifted form of the apparent hydrostatic compressive strength ($$\bar{X}/3$$) of the partially saturated
material, and $$\bar{\kappa}$$ is the branch point at which the cap function $$F_c$$ starts decreasing until it
reaches the hydrostatic strength point ($$\bar{X}$$):
<div>
$$
    \bar{\kappa} = 3\bar{p}_\text{eff}^\text{peak} - (3\bar{p}_\text{eff}^\text{peak} - \bar{X}_\text{eff}) R_c
$$
</div>
{:.notice}
where $$\bar{p}_\text{eff}^\text{peak}$$ is the maximum hydrostatic tensile stress that the material can
support and $$R_c$$ is a cap ratio.

##### The non-hardening return algorithm #####
We use the closest-point return approach discussed in the previous articles in this series to
compute a return to the yield surface while keeping all internal variables fixed.
The pseudocsode of the algorithm is listed below:
<div>
$$
  \begin{align}
    &\text{Require:} && \boldsymbol{\sigma}^k, \delta\boldsymbol{\varepsilon}, X^k, K^k, G^k,
                          (\bar{p}^w)^k, \boldsymbol{s}^\text{trial}, \sqrt{J_2^\text{trial}},
                          r^\text{trial}, z_\text{eff}^\text{trial}, \\
    &                &&  a_1, a_2, a_3, a_4, I_1^\text{peak}, R_c \\
    &\text{Procedure:} && \text{nonHardeningReturn} \\
    &   1.\quad &&   r'_\text{trial} \leftarrow r^\text{trial} \sqrt{\cfrac{3K^k}{2G^k}}
                \qquad \mbox{Transform the trial $r$ coordinate} \\
    &   2.\quad &&   X_\text{eff}^k \leftarrow X^k + 3 (\bar{p}^w)^k \\
    &   3.\quad &&   z_\text{eff}^\text{closest}, r'_\text{closest} \leftarrow
           \text{getClosestPoint}(K^k, G^k, X_\text{eff}^k, a_1, a_2, a_3, a_4,\\
    &      &&   \qquad \qquad\qquad I_1^\text{peak}, R_c, z_\text{eff}^\text{trial}, r'_\text{trial}) \\
    &   4.\quad &&   I_1^{\text{closest}} \leftarrow \sqrt{3} z_\text{eff}^\text{closest} - 3 (\bar{p}^w)^k,
              \sqrt{J_2^\text{closest}} \leftarrow \sqrt{\frac{G^k}{3K^k}}\,r'_\text{closest} \\
    &   5.\quad &&  \text{If} {\sqrt{J_2^\text{trial}} > 0} \\
    &      6. &&  \qquad \boldsymbol{\sigma}^\text{fixed} = \tfrac{1}{3} I_1^{\text{closest}} \mathbf{I} +
          \frac{\sqrt{J_2^\text{closest}}}{\sqrt{J_2^\text{trial}}}\,\boldsymbol{s}^\text{trial} \\
    &         &&  \qquad \mbox{Compute updated total stress} \\
    &   7.\quad &&  \text{Else} \\
    &      8. &&  \qquad \boldsymbol{\sigma}^\text{fixed} = \tfrac{1}{3} I_1^{\text{closest}} \mathbf{I} +
                  \boldsymbol{s}^\text{trial}\\
    &         && \qquad
                \mbox{Compute updated total stress from hydrostatic trial stress} \\
    &   9.\quad &&  \text{EndIf} \\
    &   10.\quad &&  \delta\boldsymbol{\sigma}_\text{fixed} \leftarrow \boldsymbol{\sigma}^\text{fixed} -
              \boldsymbol{\sigma}^k \qquad \mbox{Compute stress increment} \\
    &   11.\quad &&  \delta\boldsymbol{\sigma}_\text{fixed}^\text{iso} \leftarrow \tfrac{1}{3}
                \text{tr}(\delta\boldsymbol{\sigma}_\text{fixed}) \mathbf{I}, 
         \quad \delta\boldsymbol{\sigma}_\text{fixed}^\text{dev} \leftarrow
          \delta\boldsymbol{\sigma}_\text{fixed} - \delta\boldsymbol{\sigma}_\text{fixed}^\text{iso} \\
    &   12.\quad &&  \delta\boldsymbol{\varepsilon}^{p,\text{fixed}} = \delta\boldsymbol{\varepsilon} -
             \frac{1}{3K^k}\,\delta\boldsymbol{\sigma}_\text{fixed}^\text{iso} -
             \frac{1}{2G^k}\,\delta\boldsymbol{\sigma}_\text{fixed}^\text{dev} \\
    &    &&    \quad\mbox{Compute plastic strain increment} \\
    &   13.\quad &&  \text{Return} \qquad \boldsymbol{\sigma}^\text{fixed}, \delta\boldsymbol{\varepsilon}^{p,\text{fixed}} \\
  \end{align}
$$
</div>
{:.notice--info}

##### The closest point algorithm #####
Let us look at the actual implementation to get a feel for how the closest point can be found
geometrically.

{% highlight cpp %}
bool
YieldCond_Arena::getClosestPoint(const ModelState_Arena* state, // The plasticity state
                                 const double& px,              // The trial stress z, r'
                                 const double& py,
                                 double& cpx,                   // The closest point z, r'
                                 double& cpy)
{
  Point pt(px, py, 0.0);
  Point closest(0.0, 0.0, 0.0);
  getClosestPointGeometricBisect(state, pt, closest);
  cpx = closest.x();
  cpy = closest.y();
  return true;
}
{% endhighlight %}

The bisection algorithm for the closest point can be implemented as shown below.

{% highlight cpp %}
void
YieldCond_MasonSand::getClosestPointGeometricBisect(const ModelState_Arena* state,
                                                    const Point& z_r_pt,
                                                    Point& z_r_closest)
{
  // Get the particle specific internal variables from the model state
  // Store in a local struct
  d_local.PEAKI1 = state->yieldParams.at("PEAKI1");
  ................

  std::vector<double> limitParameters =
    computeModelParameters(d_local.PEAKI1, d_local.FSLOPE, d_local.STREN, d_local.YSLOPE);
  d_local.a1 = limitParameters[0];
  d_local.a2 = limitParameters[1];
  ................

  // Get the plastic internal variables from the model state
  double pbar_w = state->pbar_w;
  double X_eff = state->capX + 3.0*pbar_w;

  // Compute kappa
  double I1_diff = d_local.PEAKI1 - X_eff;
  double kappa =  d_local.PEAKI1 - d_local.CR*I1_diff;

  // Get the bulk and shear moduli and compute sqrt(3/2 K/G)
  double sqrtKG = std::sqrt(1.5*state->bulkModulus/state->shearModulus);

  // Compute diameter of yield surface in z-r space
  double sqrtJ2_diff = 2.0*evalYieldConditionMax(state);
  double yield_surf_dia_zrprime = std::max(I1_diff*one_sqrt_three, sqrtJ2_diff*sqrt_two*sqrtKG);
  double dist_to_trial_zr = std::sqrt(z_r_pt.x()*z_r_pt.x() + z_r_pt.y()*z_r_pt.y());
  double dist_dia_ratio = dist_to_trial_zr/yield_surf_dia_zrprime;

  // Set the number of points used to discretize the yield surface
  int num_points = std::max(5, (int) std::ceil(std::log(dist_dia_ratio)));

   // Set up I1 limits
  double I1eff_min = X_eff;
  double I1eff_max = d_local.PEAKI1;

  // Set up bisection
  double eta_lo = 0.0, eta_hi = 1.0;

  // Set up mid point
  double I1eff_mid = 0.5*(I1eff_min + I1eff_max);
  double eta_mid = 0.5*(eta_lo + eta_hi);

  // Do bisection
  int iters = 1;
  double TOLERANCE = 1.0e-10;
  std::vector<Point> z_r_points;
  std::vector<Point> z_r_segments;
  std::vector<Point> z_r_segment_points;
  Point z_r_closest_old;
  z_r_closest_old.x(std::numeric_limits<double>::max());
  z_r_closest_old.y(std::numeric_limits<double>::max());
  z_r_closest_old.z(0.0);

  while (std::abs(eta_hi - eta_lo) > TOLERANCE) {

    // Get the yield surface points
    z_r_points.clear();
    getYieldSurfacePointsAll_RprimeZ(X_eff, kappa, sqrtKG, I1eff_min, I1eff_max,
                                     num_points, z_r_points);

    // Find the closest point
    findClosestPoint(z_r_pt, z_r_points, z_r_closest);

    // Compute I1 for the closest point
    double I1eff_closest = sqrt_three*z_r_closest.x();

    // If (I1_closest < I1_mid)
    if (I1eff_closest < I1eff_mid) {
      I1eff_max = I1eff_mid;
      eta_hi = eta_mid;
    } else {
      I1eff_min = I1eff_mid;
      eta_lo = eta_mid;
    }

    I1eff_mid = 0.5*(I1eff_min + I1eff_max);
    eta_mid = 0.5*(eta_lo + eta_hi);

    // Distance to old closest point
    if (iters > 10 && (z_r_closest - z_r_closest_old).length2() < 1.0e-16) {
      break;
    }
    z_r_closest_old = z_r_closest;

    ++iters;
  }

  return;
}
{% endhighlight %}

The yield surface points are computed using

{% highlight cpp %}
void
YieldCond_MasonSand::getYieldSurfacePointsAll_RprimeZ(const double& X_eff,
                                                      const double& kappa,
                                                      const double& sqrtKG,
                                                      const double& I1eff_min,
                                                      const double& I1eff_max,
                                                      const int& num_points,
                                                      std::vector<Point>& z_r_vec)
{
  // Compute z_eff and r'
  computeZeff_and_RPrime(X_eff, kappa, sqrtKG, I1eff_min, I1eff_max, num_points, z_r_vec);

  return;
}
void
YieldCond_MasonSand::computeZeff_and_RPrime(const double& X_eff,
                                            const double& kappa,
                                            const double& sqrtKG,
                                            const double& I1eff_min,
                                            const double& I1eff_max,
                                            const int& num_points,
                                            std::vector<Point>& z_r_vec)
{
  // Set up points
  double rad = 0.5*(d_local.PEAKI1 - X_eff);
  double cen = 0.5*(d_local.PEAKI1 + X_eff);
  double theta_max = std::acos(std::max((I1eff_min - cen)/rad, -1.0));
  double theta_min = std::acos(std::min((I1eff_max - cen)/rad, 1.0));
  std::vector<double> theta_vec;
  linspace(theta_min, theta_max, num_points, theta_vec);

  for (auto theta : theta_vec) {
    double I1_eff = std::max(cen + rad*std::cos(theta), X_eff);


    // Compute F_f
    double Ff = d_local.a1 - d_local.a3*std::exp(d_local.a2*I1_eff) - d_local.a4*(I1_eff);
    double Ff_sq = Ff*Ff;

    // Compute Fc
    double Fc_sq = 1.0;
    if (I1_eff < kappa) {
      double ratio = (kappa - I1_eff)/(kappa - X_eff);
      Fc_sq = std::max(1.0 - ratio*ratio, 0.0);
    }

    // Compute J2
    double J2 = Ff_sq*Fc_sq;
    z_r_vec.push_back(Point(I1_eff/sqrt_three,
                            std::sqrt(2.0*J2)*sqrtKG, 0.0));
  }

  return;
}
{% endhighlight %}

And, finally, the actual geometric closest point algorithm is

{% highlight cpp %}
void
YieldCond_MasonSand::findClosestPoint(const Point& p,
                                      const std::vector<Point>& poly,
                                      Point& min_p)
{
  double TOLERANCE_MIN = 1.0e-12;
  std::vector<Point> XP;

  // Loop through the segments of the polyline
  auto iterStart = poly.begin();
  auto iterEnd   = poly.end();
  auto iterNext = iterStart;
  ++iterNext;
  for ( ; iterNext != iterEnd; ++iterStart, ++iterNext) {
    Point start = *iterStart;
    Point next  = *iterNext;

    // Find shortest distance from point to the polyline line
    Vector m = next - start;
    Vector n = p - start;
    if (m.length2() < TOLERANCE_MIN) {
      XP.push_back(start);
    } else {
      const double t0 = Dot(m, n) / Dot(m, m);
      if (t0 <= 0.0) {
        XP.push_back(start);
      } else if (t0 >= 1.0) {
        XP.push_back(next);
      } else {
        // Shortest distance is inside segment; this is the closest point
        min_p = m * t0 + start;
        XP.push_back(min_p);
      }
    }
  }

  double min_d = std::numeric_limits<double>::max();
  for (const auto& xp :  XP) {
    // Compute distance sq
    double dSq = (p - xp).length2();
    if (dSq < min_d) {
      min_d = dSq;
      min_p = xp;
    }
  }

  return;
}
{% endhighlight %}

##### An animation of the closest point algorithm #####
<div class="yield-surf-canvas">
</div>


#### Remarks ####
This geometric algorithm is remarkably accurate and avoids complications associated with
computing gradients in the transformed space.  In the next article we will discuss
how the animation in this was created.


<script src="https://d3js.org/d3.v4.min.js"></script>
<script src="{{ site.baseurl }}/assets/js/yieldsurface.js"></script>
<script>
  d3.json("{{ site.baseurl }}/assets/json/yieldSurfData.json", drawYieldSurface);
</script>

