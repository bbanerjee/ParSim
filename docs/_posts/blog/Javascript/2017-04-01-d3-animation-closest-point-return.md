---
layout: posts
title:  "Creating an animation with d3.js"
subheadline: "Biswajit Banerjee"
description: "How the closest-point return animation was created"
date:  2017-04-01 10:30:00
categories:
    - Javascript
    - D3JS
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

##### Introduction #####
In [Part 8]({{site.url }}/mechanics/plasticity/algorithm/geometric-closest-point-return/) of the
series on return algorithms for plasticity we animated the closest-point return algorithm as
shown below.

<div class="yield-surf-canvas">
</div>

An animated GIF of the plot can be generated using the following `R` script.

{% highlight r %}
require("ggplot2")
require("animation")
require("latex2exp")

#------------------------------------------------------
# Plot single iteration
#------------------------------------------------------
plotYieldSurface <- function(zrprime_data) {

  zrprime_trial_closest = zrprime_data[which(zrprime_data$Label == "trial" |
                                             zrprime_data$Label == "closest"),]
  plt = ggplot() + 
        geom_path(data = zrprime_data, 
                  aes(x = z*1.0e-6, y = rprime*1.0e-6, group=Label, color = Label),
                      size = 1)+
        geom_point(data = zrprime_data, 
                   aes(x = z*1.0e-6, y = rprime*1.0e-6, group=Label, color = Label)) +
        geom_line(data = zrprime_trial_closest,
                  aes(x = z*1.0e-6, y = rprime*1.0e-6),
                  color = "red", linetype = 1,
                  size = 1) +
        xlab(TeX(paste("$z = I_1/\\sqrt{3}$", "MPa"))) +
        ylab(TeX(paste("$r' = \\beta\\sqrt{3K/2G}\\sqrt{2J_2}$", "MPa"))) +
        coord_fixed() +
        theme_bw() + 
        theme(legend.justification=c(0,1), legend.position=c(0,1),
            plot.title = element_text(size = 14),
            axis.title.x = element_text(size=16),
            axis.title.y = element_text(size=16),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14),
            legend.text = element_text(size=12))
  print(plt)
}

#------------------------------------------------------
# Animate the plots
#------------------------------------------------------
animateIterations <- function(yieldSurfData, num_iterations, consistency_iter) {
  lapply(seq(1, num_iterations, 1), 
         function(iteration) {
           plotYieldSurface(extractIteration(yieldSurfData, iteration, consistency_iter))
         })
}

#------------------------------------------------------
# Function to create animated gif
#------------------------------------------------------
createGIFIterations <- function(yieldSurfData, consistency_iter) {

  outputGIFFile = paste0("closest_point_iter", ".gif")

  # Compute number of iterations
  num_iterations = length(unique(yieldSurfData$Iteration))

  # Save as animation
  ani.options(ani.height = 600, ani.width = 600)
  saveGIF(animateIterations(yieldSurfData, num_iterations, consistency_iter), 
          interval=1.0, 
          movie.name=outputGIFFile)
}

#-------------------------------------------------------------------------
# Compute the full yield surface and create data frame
#-------------------------------------------------------------------------
ComputeFullYieldSurface <- function(yieldParams, capX, pbar_w, K, G, num_points,
                                    z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                                    iteration, consistency_iter) {

  # Get yield parameters
  yieldParams = unlist(yieldParams)
  PEAKI1 = as.numeric(yieldParams['PEAKI1'])
  FSLOPE = as.numeric(yieldParams['FSLOPE'])
  STREN  = as.numeric(yieldParams['STREN'])
  YSLOPE = as.numeric(yieldParams['YSLOPE'])
  CR     = as.numeric(yieldParams['CR'])

  # Set up constants
  a1 = STREN
  a2 = (FSLOPE-YSLOPE)/(STREN-YSLOPE*PEAKI1)
  a3 = (STREN-YSLOPE*PEAKI1)*exp(-a2*PEAKI1)
  a4 = YSLOPE

  # Compute kappa
  X_eff = capX + 3*pbar_w
  kappa = PEAKI1 - CR*(PEAKI1 - X_eff)

  # Create an array of I1_eff values
  I1eff_min = X_eff;
  I1eff_max = PEAKI1;

  rad = 0.5*(PEAKI1 - X_eff)
  cen = 0.5*(PEAKI1 + X_eff)
  theta_max = acos(max((I1eff_min - cen)/rad, -1.0))
  theta_min = acos(min((I1eff_max - cen)/rad, 1.0))
  theta_vec = seq(from = theta_min, to = theta_max, length.out = num_points)
  I1_list = cen + rad*cos(theta_vec)
  I1_eff_list = lapply(I1_list, function(val) {max(val, X_eff)})
  sqrtJ2_list = lapply(I1_eff_list,
    function(I1_eff) {
      # Compute F_f
      Ff = a1 - a3*exp(a2*I1_eff) - a4*(I1_eff)

      # Compute Fc
      Fc_sq = 1.0
      if ((I1_eff < kappa) && (X_eff <= I1_eff)) {
        ratio = (kappa - I1_eff)/(kappa - X_eff)
        Fc_sq = 1.0 - ratio^2
      }

      # Compute sqrt(J2)
      sqrtJ2 = Ff*sqrt(Fc_sq)
      return(sqrtJ2)
    })
  zrprime_data = data.frame(z = unlist(I1_eff_list)/sqrt(3), 
                            rprime = unlist(sqrtJ2_list)*sqrt(2)*sqrt(1.5*K/G), 
                            Iteration = as.factor(iteration),
                            CIteration = as.factor(consistency_iter),
                            Label="yield")
  zrprime_data = rbind(zrprime_data,
                       data.frame(z = z_r_pt[1], rprime = z_r_pt[2], 
                                  Iteration = as.factor(iteration),
                                  CIteration = as.factor(consistency_iter),
                                  Label = "trial"))
  zrprime_data = rbind(zrprime_data,
                       data.frame(z = z_r_closest[1], rprime = z_r_closest[2], 
                                  Iteration = as.factor(iteration),
                                  CIteration = as.factor(consistency_iter),
                                  Label = "closest"))
  zrprime_data = rbind(zrprime_data,
                       data.frame(z = z_r_yield_z, rprime = z_r_yield_r, 
                                  Iteration = as.factor(iteration),
                                  CIteration = as.factor(consistency_iter),
                                  Label = "polyline"))
 
  return(zrprime_data)
}

#-------------------------------------------------------------------------
# Read data and plot for animation
#-------------------------------------------------------------------------
ReadAndPlotClosestPoint <- function() {)

  num_points = 50
  consistency_iter = 1

  iteration = 1
  K = 3.98447e+07
  G = 1.32816e+07
  X = -1.87891e+07
  pbar_w = 5.24774e+06
  yieldParams = list(BETA = 1, CR = 0.5, FSLOPE = 0.355309, PEAKI1 = 996.118, STREN = 1.76382e+07, YSLOPE = 0.353622)
  z_r_pt = c(52172.3,3.60195e+06)
  z_r_closest = c(-1.15079e+06,2.02165e+06)
  z_r_yield_z = c(575.109,-167406,-607187,-1.15079e+06,-1.59057e+06,-1.75855e+06)
  z_r_yield_r = c(2.45257e-09,310133,1.12207e+06,2.02165e+06,1.72669e+06,0)
  zr_df = 
  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, num_points,
                          z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                          iteration, consistency_iter)
  iteration = 2
  K = 3.98447e+07
  G = 1.32816e+07
  X = -1.87891e+07
  pbar_w = 5.24774e+06
  yieldParams = list(BETA = 1, CR = 0.5, FSLOPE = 0.355309, PEAKI1 = 996.118, STREN = 1.76382e+07, YSLOPE = 0.353622)
  z_r_pt = c(52172.3,3.60195e+06)
  z_r_closest = c(-1.15079e+06,2.02165e+06)
  z_r_yield_z = c(-878986,-1.15079e+06,-1.39598e+06,-1.59057e+06,-1.7155e+06,-1.75855e+06)
  z_r_yield_r = c(1.62388e+06,2.02165e+06,2.08595e+06,1.72669e+06,979052,0)
  zr_df = rbind(zr_df,
  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, num_points,
                          z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                          iteration, consistency_iter))
  iteration = 3
  K = 3.98447e+07
  G = 1.32816e+07
  X = -1.87891e+07
  pbar_w = 5.24774e+06
  yieldParams = list(BETA = 1, CR = 0.5, FSLOPE = 0.355309, PEAKI1 = 996.118, STREN = 1.76382e+07, YSLOPE = 0.353622)
  z_r_pt = c(52172.3,3.60195e+06)
  z_r_closest = c(-1.15079e+06,2.02165e+06)
  z_r_yield_z = c(-878986,-970925,-1.06186e+06,-1.15079e+06,-1.23674e+06,-1.31877e+06)
  z_r_yield_r = c(1.62388e+06,1.7838e+06,1.91864e+06,2.02165e+06,2.08688e+06,2.10948e+06)
  zr_df = rbind(zr_df,
  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, num_points,
                          z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                          iteration, consistency_iter))
  iteration = 4
  K = 3.98447e+07
  G = 1.32816e+07
  X = -1.87891e+07
  pbar_w = 5.24774e+06
  yieldParams = list(BETA = 1, CR = 0.5, FSLOPE = 0.355309, PEAKI1 = 996.118, STREN = 1.76382e+07, YSLOPE = 0.353622)
  z_r_pt = c(52172.3,3.60195e+06)
  z_r_closest = c(-1.18969e+06,2.05584e+06)
  z_r_yield_z = c(-1.09888e+06,-1.14468e+06,-1.18969e+06,-1.2338e+06,-1.27687e+06,-1.31877e+06)
  z_r_yield_r = c(1.96539e+06,2.01563e+06,2.05584e+06,2.0853e+06,2.10336e+06,2.10948e+06)
  zr_df = rbind(zr_df,
  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, num_points,
                          z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                          iteration, consistency_iter))
  iteration = 5
  K = 3.98447e+07
  G = 1.32816e+07
  X = -1.87891e+07
  pbar_w = 5.24774e+06
  yieldParams = list(BETA = 1, CR = 0.5, FSLOPE = 0.355309, PEAKI1 = 996.118, STREN = 1.76382e+07, YSLOPE = 0.353622)
  z_r_pt = c(52172.3,3.60195e+06)
  z_r_closest = c(-1.18723e+06,2.05389e+06)
  z_r_yield_z = c(-1.09888e+06,-1.12123e+06,-1.14342e+06,-1.16542e+06,-1.18723e+06,-1.20882e+06)
  z_r_yield_r = c(1.96539e+06,1.99102e+06,2.01438e+06,2.03536e+06,2.05389e+06,2.06989e+06)
  zr_df = rbind(zr_df,
  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, num_points,
                          z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                          iteration, consistency_iter))
  iteration = 6
  K = 3.98447e+07
  G = 1.32816e+07
  X = -1.87891e+07
  pbar_w = 5.24774e+06
  yieldParams = list(BETA = 1, CR = 0.5, FSLOPE = 0.355309, PEAKI1 = 996.118, STREN = 1.76382e+07, YSLOPE = 0.353622)
  z_r_pt = c(52172.3,3.60195e+06)
  z_r_closest = c(-1.18699e+06,2.05371e+06)
  z_r_yield_z = c(-1.15385e+06,-1.16495e+06,-1.176e+06,-1.18699e+06,-1.19794e+06,-1.20882e+06)
  z_r_yield_r = c(2.0246e+06,2.03493e+06,2.04464e+06,2.05371e+06,2.06213e+06,2.06989e+06)
  zr_df = rbind(zr_df,
  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, num_points,
                          z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                          iteration, consistency_iter))
  iteration = 7
  K = 3.98447e+07
  G = 1.32816e+07
  X = -1.87891e+07
  pbar_w = 5.24774e+06
  yieldParams = list(BETA = 1, CR = 0.5, FSLOPE = 0.355309, PEAKI1 = 996.118, STREN = 1.76382e+07, YSLOPE = 0.353622)
  z_r_pt = c(52172.3,3.60195e+06)
  z_r_closest = c(-1.18686e+06,2.0536e+06)
  z_r_yield_z = c(-1.18133e+06,-1.18686e+06,-1.19237e+06,-1.19787e+06,-1.20335e+06,-1.20882e+06)
  z_r_yield_r = c(2.04911e+06,2.0536e+06,2.05792e+06,2.06208e+06,2.06607e+06,2.06989e+06)
  zr_df = rbind(zr_df,
  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, num_points,
                          z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                          iteration, consistency_iter))
  iteration = 8
  K = 3.98447e+07
  G = 1.32816e+07
  X = -1.87891e+07
  pbar_w = 5.24774e+06
  yieldParams = list(BETA = 1, CR = 0.5, FSLOPE = 0.355309, PEAKI1 = 996.118, STREN = 1.76382e+07, YSLOPE = 0.353622)
  z_r_pt = c(52172.3,3.60195e+06)
  z_r_closest = c(-1.18684e+06,2.05359e+06)
  z_r_yield_z = c(-1.18133e+06,-1.18409e+06,-1.18684e+06,-1.18959e+06,-1.19234e+06,-1.19508e+06)
  z_r_yield_r = c(2.04911e+06,2.05137e+06,2.05359e+06,2.05576e+06,2.05789e+06,2.05999e+06)
  zr_df = rbind(zr_df,
  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, num_points,
                          z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                          iteration, consistency_iter))
  iteration = 9
  K = 3.98447e+07
  G = 1.32816e+07
  X = -1.87891e+07
  pbar_w = 5.24774e+06
  yieldParams = list(BETA = 1, CR = 0.5, FSLOPE = 0.355309, PEAKI1 = 996.118, STREN = 1.76382e+07, YSLOPE = 0.353622)
  z_r_pt = c(52172.3,3.60195e+06)
  z_r_closest = c(-1.18683e+06,2.05358e+06)
  z_r_yield_z = c(-1.18133e+06,-1.18271e+06,-1.18409e+06,-1.18546e+06,-1.18683e+06,-1.18821e+06)
  z_r_yield_r = c(2.04911e+06,2.05025e+06,2.05137e+06,2.05248e+06,2.05358e+06,2.05467e+06)
  zr_df = rbind(zr_df,
  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, num_points,
                          z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                          iteration, consistency_iter))
  iteration = 10
  K = 3.98447e+07
  G = 1.32816e+07
  X = -1.87891e+07
  pbar_w = 5.24774e+06
  yieldParams = list(BETA = 1, CR = 0.5, FSLOPE = 0.355309, PEAKI1 = 996.118, STREN = 1.76382e+07, YSLOPE = 0.353622)
  z_r_pt = c(52172.3,3.60195e+06)
  z_r_closest = c(-1.18658e+06,2.05338e+06)
  z_r_yield_z = c(-1.18477e+06,-1.18546e+06,-1.18615e+06,-1.18683e+06,-1.18752e+06,-1.18821e+06)
  z_r_yield_r = c(2.05192e+06,2.05248e+06,2.05303e+06,2.05358e+06,2.05412e+06,2.05467e+06)
  zr_df = rbind(zr_df,
  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, num_points,
                          z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                          iteration, consistency_iter))
  iteration = 11
  K = 3.98447e+07
  G = 1.32816e+07
  X = -1.87891e+07
  pbar_w = 5.24774e+06
  yieldParams = list(BETA = 1, CR = 0.5, FSLOPE = 0.355309, PEAKI1 = 996.118, STREN = 1.76382e+07, YSLOPE = 0.353622)
  z_r_pt = c(52172.3,3.60195e+06)
  z_r_closest = c(-1.18649e+06,2.0533e+06)
  z_r_yield_z = c(-1.18649e+06,-1.18683e+06,-1.18718e+06,-1.18752e+06,-1.18786e+06,-1.18821e+06)
  z_r_yield_r = c(2.0533e+06,2.05358e+06,2.05385e+06,2.05412e+06,2.0544e+06,2.05467e+06)
  zr_df = rbind(zr_df,
  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, num_points,
                          z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                          iteration, consistency_iter))
  iteration = 12
  K = 3.98447e+07
  G = 1.32816e+07
  X = -1.87891e+07
  pbar_w = 5.24774e+06
  yieldParams = list(BETA = 1, CR = 0.5, FSLOPE = 0.355309, PEAKI1 = 996.118, STREN = 1.76382e+07, YSLOPE = 0.353622)
  z_r_pt = c(52172.3,3.60195e+06)
  z_r_closest = c(-1.18649e+06,2.0533e+06)
  z_r_yield_z = c(-1.18649e+06,-1.18666e+06,-1.18683e+06,-1.187e+06,-1.18718e+06,-1.18735e+06)
  z_r_yield_r = c(2.0533e+06,2.05344e+06,2.05358e+06,2.05371e+06,2.05385e+06,2.05399e+06)
  zr_df = rbind(zr_df,
  ComputeFullYieldSurface(yieldParams, X, pbar_w, K, G, num_points,
                          z_r_pt, z_r_closest, z_r_yield_z, z_r_yield_r,
                          iteration, consistency_iter))
  createGIFIterations(zr_df, consistency_iter)
}

#-------------------------------------------------------------------------
# Call function to read data and plot
#-------------------------------------------------------------------------
ReadAndPlotClosestPoint()

{% endhighlight %}

Let us explore how the animation at the top of the docs was produced using [d3.js](https://d3js.org/).

#### Input data ####
The input data that can be seen in the `R` function `ReadAndPlotClosestPoint` were generated
during a simulation using print statements.  These have to be converted into a form that can
be read easily by Javascript.  For our simulation we converted these into JSON and created the
file `yieldSurfData.json` that contains:

{% highlight json %}
{
  "num_points":50,
  "consistency_iter":1,
  "data":[
    {
      "iteration":1,
      "K":3.98447e+07,
      "G":1.32816e+07,
      "X":-1.87891e+07,
      "pbar_w":5.24774e+06,
      "BETA":1,
      "CR":0.5,
      "FSLOPE":0.355309,
      "PEAKI1":996.118,
      "STREN":1.76382e+07,
      "YSLOPE":0.353622,
      "z_r_pt":[
        52172.3,
        3.60195e+06
      ],
      "z_r_closest":[
        -1.15079e+06,
        2.02165e+06
      ],
      "z_r_yield_z":[
        575.109,
        -167406,
        -607187,
        -1.15079e+06,
        -1.59057e+06,
        -1.75855e+06
      ],
      "z_r_yield_r":[
        2.45257e-09,
        310133,
        1.12207e+06,
        2.02165e+06,
        1.72669e+06,
        0
      ]
    },
    .....
  ]
}
{% endhighlight %}

#### The HTML file ####
In the HTML file, I added the following to allow the animation to be added to the DOM.

{% highlight html %}
  <body>
     ...
     <!-- Create div where the svg canvas will be added -->
     <div class="yield-surf-canvas"> </div>
     <!-- Load D3.js -->
     <script src="https://d3js.org/d3.v4.min.js"></script>
     <!-- Load javascript for animating yield surface -->
     <script src="{{ site.url }}/assets/js/yieldsurface.js"></script>
     <script>
       // Read JSON file and then call javascript function to animate yield surface 
       d3.json("{{ site.url }}/assets/json/yieldSurfData.json", drawYieldSurface);
     </script>
  </body>
{% endhighlight %}

#### The Javascript code ####
Let us now explore the Javascript code used to generate the animation.

##### Generating points on the yield surface ####
We first generate points on the yield surface in $$z-r'$$-space.  These are displayed in blue in the figure.
The function `computeYieldSurface` takes an input the data for a single iteration (read from the
JSON file) and the number of points to be used to discretize the yield surface.

{% highlight js %}
function computeYieldSurface(iterData, numPts) {

  // Get the moduli
  let K = iterData['K'];
  let G = iterData['G'];

  // Get yield parameters
  let PEAKI1 = iterData['PEAKI1'];
  let FSLOPE = iterData['FSLOPE'];
  let STREN  = iterData['STREN'];
  let YSLOPE = iterData['YSLOPE'];
  let CR     = iterData['CR'];

  // Set up constants
  let a1 = STREN;
  let a2 = (FSLOPE-YSLOPE)/(STREN-YSLOPE*PEAKI1);
  let a3 = (STREN-YSLOPE*PEAKI1)*Math.exp(-a2*PEAKI1);
  let a4 = YSLOPE;

  // Compute kappa
  let X_eff = iterData['X'] + 3.0*iterData['pbar_w'];
  let kappa = PEAKI1 - CR*(PEAKI1 - X_eff);

  // Create an array of I1_eff values
  let I1eff_min = X_eff;
  let I1eff_max = PEAKI1;

  let rad = 0.5*(PEAKI1 - X_eff);
  let cen = 0.5*(PEAKI1 + X_eff);
  let theta_max = Math.acos(Math.max((I1eff_min - cen)/rad, -1.0));
  let theta_min = Math.acos(Math.min((I1eff_max - cen)/rad, 1.0));
  let theta_inc = (theta_max-theta_min)/numPts;
  let theta_vec = d3.range(theta_min, theta_max+theta_inc, theta_inc);
  let I1_list = theta_vec.map(function(theta) {return cen + rad*Math.cos(theta);});
  let I1_eff_list = I1_list.map(function(I1) {return Math.max(I1, X_eff);});

  // Create an array of sqrtJ2 values
  let sqrtJ2_list = I1_eff_list.map(
          function(I1_eff) {
            // Compute F_f
            let Ff = a1 - a3*Math.exp(a2*I1_eff) - a4*(I1_eff);

            // Compute Fc
            let Fc_sq = 1.0;
            if ((I1_eff < kappa) && (X_eff <= I1_eff)) {
              let ratio = (kappa - I1_eff)/(kappa - X_eff);
              Fc_sq = 1.0 - ratio*ratio;
            }

            // Compute sqrt(J2)
            let sqrtJ2 = Ff*Math.sqrt(Fc_sq);
            return(sqrtJ2);
          });

  // Compute z and r'
  let z_list = I1_eff_list.map(function(I1_eff) {return I1_eff/Math.sqrt(3.0);});
  let rprime_list = sqrtJ2_list.map(function(sqrtJ2) {return sqrtJ2*Math.sqrt(3.0*K/G);});

  // Create zipped array
  let zrprime_data = d3.zip(z_list, rprime_list);
  return zrprime_data;
}
{% endhighlight %}

The procedure is identical to that used in the `R` script, except for two differences:

* the `d3.range` function is used to create the vector of $$\theta$$ values.  This function
  differs from the `R` function `seq` in that the end point of the range is not included in
  the sequence that is produced.  Therefore, the set of points produced by `d3.js` is not
  identical to that produced by `R`.
* the function returns a zipped array created with `d3.zip` instead of a `R` dataframe.

##### Drawing the yield surface and closest-point projection ####
We use the `drawYieldSurface` function to plot and animate the yield surface and the closest-point
projection.  Let us look at the full code first and then explore some of the details.

{% highlight js %}
function drawYieldSurface(data) {

  // Create svg
  var svgsize = {x: 600, y: 600};
  var margin = {top: 50, right: 50, bottom: 50, left: 50};
  var width = svgsize.x - margin.left - margin.right;
  var height = svgsize.y - margin.top - margin.bottom;
  var svg = d3.select(".yield-surf-canvas")
              .append("svg")
                 .attr("width", svgsize.x)
                 .attr("height", svgsize.y);

  // Create chart area
  var chart = svg.append("g")
                    .attr("class", "chart")
                    .attr("transform", 
                          "translate(" + margin.left + "," + margin.top + ")");

  // Process the input data
  let numPts = data['num_points'];
  let numIter = data.data.length;

  // Compute yield surface in z and r' from first element of data
  // (Note: the yield surface does not change at each iteration)
  let iterData = data.data[0];
  let zrprime = computeYieldSurface(iterData, numPts);

  // Get the coordinates of the trial stress
  let ztrial = iterData.z_r_pt[0];
  let rprimetrial = iterData.z_r_pt[1];

  // Compute min max of domain
  let zmin = d3.min(zrprime, function(d) {return d[0];})
  let zmax = d3.max(zrprime, function(d) {return d[0];})
  let rprimemin = d3.min(zrprime, function(d) {return d[1];})
  let rprimemax = d3.max(zrprime, function(d) {return d[1];})
  zmin = Math.min(zmin, ztrial);
  zmax = Math.max(zmax, ztrial);
  rprimemin = Math.min(rprimemin, rprimetrial);
  rprimemax = Math.max(rprimemax, rprimetrial);

  // Create scaling functions
  var xscale = d3.scaleLinear().domain([zmin, zmax]).rangeRound([0, width]);
  var yscale = d3.scaleLinear().domain([rprimemin, rprimemax]).rangeRound([height, 0]);
  var xscaleAxis = d3.scaleLinear().domain([zmin*1.0e-6, zmax*1.0e-6]).rangeRound([0, width]);
  var yscaleAxis = d3.scaleLinear().domain([rprimemin*1.0e-6, rprimemax*1.0e-6]).rangeRound([height, 0]);

  // Create polyline generation function
  var line = d3.line()
                 .x(function(d) { return xscale(d[0]); })
                 .y(function(d) { return yscale(d[1]); });

  // Create bottom axis
  chart.append("g")
          .attr("transform", "translate(0," + height + ")")
        .call(d3.axisBottom(xscaleAxis))
        .append("text")
           .attr("fill", "#000")
           .attr("transform",
                 "translate(" + (width/2) + " ," + 
                                (margin.top - 20) + ")")
           .attr("text-anchor", "end")
           .text("z (MPa)"); 

  // Create left axis
  chart.append("g")
       .call(d3.axisLeft(yscaleAxis))
       .append("text")
          .attr("fill", "#000")
          .attr("transform", "rotate(-90)")
          .attr("y", 6)
          .attr("dy", "0.71em")
          .attr("text-anchor", "end")
          .text("r' (MPa)"); 

  // Plot yield surface
  chart.append("path")
          .attr("class", "yield_surf")
          .datum(zrprime)
          .attr("fill", "none")    
          .attr("stroke", "steelblue")
          .attr("stroke-linejoin", "round")
          .attr("stroke-linecap", "round")
          .attr("stroke-width", 1.5)
          .attr("d", line);

  // Loop through the iterations
  for (let iter = 0; iter < numIter; iter++) {

    iterData = data.data[iter];

    // Create group for each iteration
    let iterGroup = chart.append("g")
                            .attr("class", "iteration"+iter)
                            .attr("opacity", 0);

    // Create group for the yield surface points + circles
    let yieldSurfGroup = iterGroup.append("g") 
                                     .attr("class", "approx_yield_surf"+iter);

    // Get the coordinates of polyline approximation
    let zpoly = iterData.z_r_yield_z;
    let rprimepoly = iterData.z_r_yield_r;
    let zrprimepoly = d3.zip(zpoly, rprimepoly);

    // Add the yield surface polyline to the svg
    yieldSurfGroup.append("path")
                  .datum(zrprimepoly)
                     .attr("fill", "none")    
                     .attr("stroke", "red")
                     .attr("stroke-linejoin", "round")
                     .attr("stroke-linecap", "round")
                     .attr("stroke-width", 1.5)
                     .attr("d", line);

    // Add circles to the yield surface polyline
    yieldSurfGroup.selectAll(".approx_yield_surf_circle"+iter)
                  .data(zrprimepoly)
                  .enter()
                  .append("circle")
                     .attr("r", 2)
                     .attr("cx", function(d) {return xscale(d[0]);})
                     .attr("cy", function(d) {return yscale(d[1]);})
                     .attr("fill", "blue")
                     .attr("stroke", "black");

    // Create group for the projection line and points
    let projLineGroup = iterGroup.append("g") 
                                    .attr("class", "projection_line"+iter);

    // Get the coordinates of the trial stress
    let ztrial = iterData.z_r_pt[0];
    let rprimetrial = iterData.z_r_pt[1];

    // Get the coordinates of the closest point stress
    let zclosest = iterData.z_r_closest[0];
    let rprimeclosest = iterData.z_r_closest[1];

    // Set up projection line
    let zrprimeproj = [[ztrial, rprimetrial],[zclosest, rprimeclosest]];

    // Add the projection line to the svg
    projLineGroup.append("path")
                 .datum(zrprimeproj)
                   .attr("fill", "none")    
                   .attr("stroke", "green")
                   .attr("stroke-linejoin", "round")
                   .attr("stroke-linecap", "round")
                   .attr("stroke-width", 1.5)
                   .attr("d", line);

    // Add the triangles at the end points of the projection line
    projLineGroup.selectAll(".projection_line_point"+iter)
                 .data(zrprimeproj)
                 .enter()
                 .append("path")
                    .attr("class", "projection_line_point"+iter)
                    .attr("d", d3.symbol().type(d3.symbolTriangle))
                    .attr("transform", 
                      function(d) {return "translate(" + xscale(d[0]) + "," + yscale(d[1]) + ")";});

  } // End for loop over iterations

  // Function to do the animation
  function drawAnimation() {
    var iterationID = 0;

    // Function to draw the surface for the next iteration
    function drawNextSurf() {

      d3.active(this)
          .attr("opacity", 1)
        .transition()
          .duration(1000)
          .attr("opacity", 0);
      iterationID++;
      if (iterationID === numIter) {
        iterationID = 0;
      }
      chart.select(".iteration"+iterationID)
           .transition()
             .delay(1000)
             .duration(1000)
           .on("start", drawNextSurf);            

    } // end function drawNextSurf

    // Start the animation
    chart.select(".iteration0")
         .transition()
           .duration(1000)
         .on("start", drawNextSurf);

  } // end function drawAnimation

  // Do the animation
  drawAnimation();

} 
{% endhighlight %}

Let us now look at some of the details of the code.

###### Creating the SVG ######
First we select the `<div>` in
our HTML file using `d3.select` and add a SVG canvas to it:
{% highlight js %}
var svg = d3.select(".yield-surf-canvas")
            .append("svg")
              .attr("width", svgsize.x)
              .attr("height", svgsize.y);
{% endhighlight %}

###### Creating the canvas group ######
Then we create a subset of the SVG canvas as the region where the plot will be
produced and identify it using a group `g`.  This area has to be translated
so that the top left hand corner is at the right position.
{% highlight js %}
var chart = svg.append("g")
                 .attr("class", "chart")
                 .attr("transform", 
                       "translate(" + margin.left + "," + margin.top + ")");
{% endhighlight %}

###### Creating map from real to canvas coordinates ######
We then compute the true yield surface and find the extents of the domain of
the plot based on the yield surface and the trial stress state.  We use
these extents to create maps from real coordinates to svg coordinates.
Separate maps are also created so that the axis labels can be converted into MPa.
{% highlight js %}
var xscale = d3.scaleLinear()
                 .domain([zmin, zmax])
                 .rangeRound([0, width]);
var yscale = d3.scaleLinear()
                 .domain([rprimemin, rprimemax])
                 .rangeRound([height, 0]);
var xscaleAxis = d3.scaleLinear()
                     .domain([zmin*1.0e-6, zmax*1.0e-6])
                     .rangeRound([0, width]);
var yscaleAxis = d3.scaleLinear()
                     .domain([rprimemin*1.0e-6, rprimemax*1.0e-6])
                     .rangeRound([height, 0]);
{% endhighlight %}

###### Creating SVG polyline generator ######
To plot the polylines that represent the yield surface, its approximation, and
the closest point projection, we need a function that can convert a set of points
to the SVG representation of a polyline.  We generate that function using
`d3.line`:
{% highlight js %}
var line = d3.line()
               .x(function(d) { return xscale(d[0]); })
               .y(function(d) { return yscale(d[1]); });
{% endhighlight %}

###### Creating the axes and the yield surface ######
Next, we create the axes.  The generation of axes using `d3.js` is straightforward,
but we have to keep in mind that a translation may be required depending on the
position of the axis.

After the axes have been created, we add the polyline representing the yield surface
to the SVG.  This is done using the zipped data returned from the `computeYieldSurface`
function.  Because the yield surface does not change during closest-point search iterations,
we use the data from the first iteration to plot this curve.
{% highlight js %}
chart.append("path")
       .attr("class", "yield_surf")
     .datum(zrprime)
       .attr("fill", "none")    
       .attr("stroke", "steelblue")
       .attr("stroke-width", 1.5)
       .attr("d", line);
{% endhighlight %}
Notice that we use the `.datum` method in this case and that the attribute `d` is set
using the `line` function we had defined earlier.

###### Creating groups for each iteration ######
Now that the fixed content of the plot has been created, we are ready to create the
content that changes between iterations.  To do that, we loop through the iterations
and add a group for each iteration:
{% highlight js %}
let iterGroup = chart.append("g")
                       .attr("class", "iteration"+iter)
                       .attr("opacity", 0);
{% endhighlight %}
We set the opacity of the group to 0 so that although all the objects in the
plot have been created, none are visible when the animation starts.

###### Adding circles to the approximate surface ######
Now we are ready to add in the polyline approximation to the yield surface
in red.  The procedure is identical to the plot of the full yield surface.
We also add circles to each vertex of the polyline using
{% highlight js %}
yieldSurfGroup.selectAll(".approx_yield_surf_circle"+iter)
                .data(zrprimepoly)
                .enter()
                .append("circle")
                  .attr("class", "approx_yield_surf_circle"+iter)
                  .attr("r", 2)
                  .attr("cx", function(d) {return xscale(d[0]);})
                  .attr("cy", function(d) {return yscale(d[1]);})
                  .attr("fill", "blue")
                  .attr("stroke", "black");
{% endhighlight %}

###### Adding triangles to the projection line ######
Next we add the line that represents the closest-point projection (in green)
and add triangles to the end points for clarity.  The triangles are
added using
{% highlight js %}
projLineGroup.selectAll(".projection_line_point"+iter)
               .data(zrprimeproj)
               .enter()
               .append("path")
                 .attr("class", "projection_line_point"+iter)
                 .attr("d", d3.symbol().type(d3.symbolTriangle))
                 .attr("transform", 
                      function(d) {return "translate(" + xscale(d[0]) + "," + yscale(d[1]) + ")";});
{% endhighlight %}

###### Setting up the animation ######
Now we can terminate the loop of the iterations and proceed to setting up the animation.
We start the animation by selecting the group that has class `iteration0` and use a transition
that takes 1 second to move from one frame to the next.  We call the `drawNextSurf` method when
the animation starts.
{% highlight js %}
chart.select(".iteration0")
     .transition()
       .duration(1000)
     .on("start", drawNextSurf);
{% endhighlight %}

The `drawNextSurf` method changes the opacity of the active group to 1, waits for 1 second,
and then changes back the opacity to 0.
{% highlight js %}
d3.active(this)
    .attr("opacity", 1)
  .transition()
    .duration(1000)
    .attr("opacity", 0);
{% endhighlight %}
After that, it increments to iteration count, selects the next group, and recurses after a short delay.
{% highlight js %}
chart.select(".iteration"+iterationID)
     .transition()
       .delay(1000)
       .duration(1000)
     .on("start", drawNextSurf);
{% endhighlight %}
To allow for the animation to loop, we set the iteration index to 0 when the total number of
iterations in the data is reached.


#### Remarks ####
This `d3.js` animation took me 2 days to complete, starting from absolutely no knowledge of the
library, and with significant breaks in the workflow.  I'd suggest that, since a novice can learn and
use the library so quickly, `d3.js` should be an essential tool in the computational engineers toolbox
to make data available in the form of interactive plots instead of the static (and hard to grasp) plots
that are the norm in engineering literature.

In the next article, we will go back to the work on XML particle data that we had started earlier
and discuss ways to reading XML files using C++ code.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

<script src="https://d3js.org/d3.v4.min.js"></script>
<script src="{{ site.url }}/assets/js/yieldsurface.js"></script>
<script>
  d3.json("{{ site.url }}/assets/json/yieldSurfData.json", drawYieldSurface);
</script>

