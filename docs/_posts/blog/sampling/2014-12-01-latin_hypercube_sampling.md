---
title:  "Sampling large data sets (Part 1)"
subheadline: "Latin hypercubes"
description: "Recently, we ran into a data set that contained the positions and time-series data for a set of 50,000 sensors arranged in a cylindrical array." 
date:  2014-12-01 00:00:00
categories:
    - sampling
excerpt_separator: <!--more-->
toc: true
toc_label: "Contents"
toc_sticky: true
---

<img style="float:right;width:300px" alt="SensorArray" src="{{site.baseurl}}/assets/blogimg/SensorArray.png"/>  
The amount of data is large.  Before we analyze the whole set of data, we can find trends  
by examining a smaller set of sensors.  A smaller set is also useful if we want to
train a model. 
<!--more-->

One of extracting a smaller set from the data is to generate a [Latin hypercube sample](http://en.wikipedia.org/wiki/Latin_hypercube_sampling) from a uniform three-dimensional distribution.
The R package `lhs` can be used to generate a sample:

~~~ R
  # Install and load the lhs library
  install.packages("lhs")
  library("lhs")

  # Create a Latin Hypercube sample of coordinates between 0 and 1
  # 100 samples, 3 dimensions
  sampleCoords <- randomLHS(100, 3)
~~~

<img style="float:right;width:300px" alt="LHS samples" src="{{site.baseurl}}/assets/blogimg/SampleCoords.png"/>  
The sampling algorithm produces a set of points whose coordinates are in $$[0,1]$$.  The 
adjacent figure shows some of the points.  Blue circles indicate the distance of the point
from the plane of the screen.

To get a set of sample points that samples the sensor positions equally well, we have to
map the cube into a cylinder.  

A [simple map](http://mathproofs.blogspot.co.nz/2005/07/mapping-square-to-circle.html) that takes lines in the $$[-1,1]$$ square to ellipses in the unit circle is

$$
  r = \sqrt{x^2 + y^2 - x^2\,y^2}  \,,
$$
$$
  \theta = \tan^{-1}\left(\cfrac{y\sqrt{1-x^2/2}}{x\sqrt{1-y^2/2}}\right) \,,
$$
$$
  z = z \,.
$$

A slightly better map is the [Shirley-Chew idea](https://mediatech.aalto.fi/~jaakko/T111-5310/K2013/JGT-97.pdf) that attempts to reduce distortion during mapping.
Polar plots of the mapped LHS points show how the two maps operate.

|<img style="float:right;width:300px" alt="PolarPlots" src="{{site.baseurl}}/assets/blogimg/PolarPrime.png"/>|<img style="float:right;width:300px" alt="SensorArray" src="{{site.baseurl}}/assets/blogimg/PolarShirleyChew.png"/>|

<img style="float:right;width:300px" alt="Sensors and samples" src="{{site.baseurl}}/assets/blogimg/SensorsAndSamples.png"/>  
The sample coordinates generated in the unit disk can then be transformed so that 
they lie in the first quadrant.

Now that we have the samples, the next step is to find the sensors that are
closest to the sample points in some manner.  The next post will discuss how
that is done.

