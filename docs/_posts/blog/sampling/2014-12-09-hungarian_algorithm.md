---
title:  "Sampling large data sets (Part 2)"
subheadline: "Hungarian algorithm"
description: "In the previous post, we discussed a technique for choosing a Latin Hybercube sample within a cylindrical domain. We now need to match the sample to sensors arranged in a cylindrical array. One way of approaching this problem is to use the Hungarian algorithm. " 
date:  2014-12-09 00:00:00
categories:
    - sampling
excerpt_separator: <!--more-->
toc: true
toc_label: "Contents"
toc_sticky: true
---

<img style="float:right;width:300px" alt="Sensors and samples" src="{{site.baseurl}}/assets/blogimg/SensorsAndSamples.png"/>  
Recall that we had found a set of points that samples a cylindrical domain uniformly.  
However, these points do not match the locations of our sensors and we will have to
find the nearest sensor to each sample point.
<!--more-->

This problem is a version of the [Assignment Problem](http://en.wikipedia.org/wiki/Assignment_problem) and an elegant way of proceeding is to use the [Hungarian Algorithm](http://en.wikipedia.org/wiki/Hungarian_algorithm).  A clear description of the problem and the algorithm 
that is used to solve it can be found at [http://www.math.harvard.edu/](http://www.math.harvard.edu/archive/20_spring_05/handouts/assignment_overheads.pdf).

We have 50,000 sensors and 100 sample points.  If we want to solve the assignment problem 
directly on our desktop, the memory requirement for the complete problem is around 23 Gb
- more than the free memory available in most desktops.  So we have to simplify the 
problem to suit the resources that we have available.  

The R package `fields` has a method `rdist` for computing the Euclidean distance matrix.  We 
first find the pairwise distances between the sensors and the samples

~~~ R
  # Install and load the fields library
  install.packages("fields")
  library("fields")

  # Get the Euclidean distance matrix 
  # between the sensors and the samples
  # 50,000 sensors and 100 samples
  distances <- rdist(inputData, sampleData)

~~~

<img style="float:right;width:300px" alt="Sensors and samples" src="{{site.baseurl}}/assets/blogimg/SparseSensorsAndSamples.png"/>  
Next we define a radius of support around each sample location and find the sensors 
in that region.

~~~ R
  # Find a sparser matrix based on a radius of support
  supportRadius <- max(rdist(inputData[1:1],inputData[4:4]))
  distanceFlags <- ifelse(distances < supportRadius, 0, 1)
  sparseIndices <- which(distanceFlags == 0, arr.ind=T)
  sparseInputData <- inputData[sparseIndices[,1]]

~~~
\\
This procedure gives us the sensor locations (in red) that we care about.  Now we can 
recompute a much smaller distance matrix and apply the Hungarian algorithm to find
which sensors are closest to our sample points.

<img style="float:right;width:300px" alt="Sensors and samples" src="{{site.baseurl}}/assets/blogimg/SampleSensors.png"/>  
The `clue` package in R provides us with the required tools to solve the
<em>Linear Sum Assignment Problem</em> and the computation is 
direct:

~~~ R
# Find the distances between the input data and the sample data
distancesSparse <- rdist(sparseInputData, sampleData)

# Use Hungarian algorithm to minimize pairwise Euclidean norm
sol <- solve_LSAP(t(distancesSparse))

~~~
\\
The sample sensors can now be identified from the output of the assignment algorithm and
use in our training exercise.

If the training set is a large fraction of the total number of sensors (typically 80%-90%),
the approach that we have used becomes inefficient and other sampling techniques 
may be preferred.  We will talk about some of of these in the next part of this series.

