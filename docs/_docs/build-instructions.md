---
title: Building Vaango
layout: single
author_profile: true
permalink: "docs/build-instructions/"
sidebar:
  nav: "parsim_docs"
---
---

Instructions for building Vaango on Ubuntu are given below.

* Contents
{:toc}


## Prerequisites

### Cmake:
* You will probably need to install `Cmake` to control the software compilation process. To do that:

{% highlight sh %}
     sudo apt-get install cmake
{% endhighlight %}

### Compilers:
* You also need to have `gfortran`:

{% highlight sh %}
     sudo apt-get install gfortran
{% endhighlight %}

* and a C Compiler known as `gcc`:

{% highlight sh %}
     sudo apt-get install gcc
{% endhighlight %}

* You will need `boost` for parts of the code (and also the unit tests) to compile.

{% highlight sh %}
     sudo apt-get install libboost-all-dev
{% endhighlight %}

### MPI and XML libraries:
* the `OpenMPI` libraries:

{% highlight sh %}
     sudo apt-get install mpi-default-dev
{% endhighlight %}

* and `libxml2`:

{% highlight sh %}
     sudo apt-get install libxml2
     sudo apt-get install libxml2-dev
{% endhighlight %}

### Other libraries:
* You will also need to install the development version of zlib.

{% highlight sh %}
     sudo apt-get install zlib1g zlib1g-dev
{% endhighlight %}

* The Peridynamics code uses parts of the `Eigen3` library. You will have to install this library if you don't have it in your system:

{% highlight sh %}
     sudo apt-get install libeigen3-dev
{% endhighlight %}

* In order for related code like MPM3D_xx etc. to work, you will also need the `VTK` libraries. Use

{% highlight sh %}
     sudo apt-get install libvtk5-dev
     sudo apt-get install python-vtk tcl-vtk libvtk-java libvtk5-qt4-dev
{% endhighlight %}

* In some older versions of `Vaango`, `libpcl` is used read point data file. In order to get PCL(The Point Cloud Library) you need to run these three commands:

{% highlight sh %}
    sudo add-apt-repository ppa:v-launchpad-jochen-sprickerhof-de/pcl
    sudo apt-get update
    sudo apt-get install libpcl-all
{% endhighlight %}

## Building the executables

### The Vaango repository
After you get the code from GitHub, follow these steps:
* Go to the `Vaango` directory:

{% highlight sh %}
     cd ParSim/Vaango
{% endhighlight %}

*  The source code is in the directory `src`.

### Out-of-source optimized build
* For the optimized build, create a new directory called `opt` under `Vaango`:

{% highlight sh %}
    mkdir opt
{% endhighlight %}

* followed by:

{% highlight sh %}
    cd opt
{% endhighlight %}

* To create the make files do:

{% highlight sh %}
    cmake ../src
{% endhighlight %}

### Out-of-source debug build
* For the debug build, create a directory `dbg` under `Vaango`:

{% highlight sh %}
    mkdir dbg
    cd dbg
{% endhighlight %}

* To create the makefiles for the debug build, use

{% highlight sh %}
    cmake -DCMAKE_BUILD_TYPE=Debug ../src
{% endhighlight %}

### Unit tests
If you want the units tests to be compiled, use the alternative command

{% highlight sh %}
    cmake ../src -DBUILD_UNITS_TESTS=1
{% endhighlight %}

### Clang compiler:
If you want to used the `clang` compiler instead of `gcc`:

{% highlight sh %}
    cmake ../src -DUSE_CLANG=1
{% endhighlight %}

### Visit build:
Older versions of `Visit` required the following extra step if you want to make sure that Visit is able to read Uintah output format files (also called UDA files) you will need to use

{% highlight sh %}
   cmake ../src -DVISIT_DIR=/path/to/Visit
{% endhighlight %}

## Compiling the code:
* Next you need to compile the needed files from `src`. So enter your `opt` directory and then type : 

{% highlight sh %}
...:~/ParSim/Vaango/opt$ make -j4
{% endhighlight %}

* After this compilation you have all the executable files in your `opt` directory.
* The same process can also be used for the debug build.


