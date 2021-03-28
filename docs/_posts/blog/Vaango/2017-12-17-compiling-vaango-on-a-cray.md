---
title:  "Compiling and running the MPM code Vaango on a Cray"
subheadline: "Biswajit Banerjee"
description: "Copper and Excalibur"
date:  2017-12-16 00:30:00
categories:
    - Vaango
excerpt_separator: <!--more-->
toc: true
toc_label: "Contents"
toc_sticky: true
---


#### Introduction ####

One of the reasons I switched to `cmake` for my builds was the need to compile my
`Vaango` code on a BlueGene/Q system.  The code was previously configured using `autoconf` and
`m4`,  and because I hadn't done that implementation myself it was nightmare to change things
when I wanted to add some new feature or library.  After moving to `cmake` things have
become decidedly easier.  In this article I will discuss how I ported my code to a couple of
Cray XE6m machines and how I run the code on that machine. You can download the latest version of the code
from [GitHub](https://github.com/bbanerjee/ParSim/tree/master/Vaango).
<!--more-->

#### Authentication ####

The authentication system for the Cray machine uses [Kerberos](https://web.mit.edu/kerberos/)
and requires a [Yubikey](https://www.yubico.com/).  The typical process for authentication
is as follows:

{% highlight bash %}
1) kshell
2) kinit <username>@<domain.name>
2a) <password>
2b) authenticate with Yubikey
3) klist
4) kssh <username>@<machine-name>.<domain-name>
{% endhighlight %}

Here `kssh` is a Kerberized version of `ssh`.

#### Building the `Vaango` code on `Copper` (Cray XE6m) ####

The `Vaango` code is written in C++ and has over the past two years slowly added
several C++11 and more recently C++14 features.  I first had to make sure that
`wget` and `git` were available on the Cray machines.

##### Downloading the code #####

To download the code I used the standard procedure
{% highlight bash %}
git clone --recursive https://github.com/bbanerjee/ParSim.git
{% endhighlight %}

I then made sure that my `googletest` and `json` submodules had been downloaded
correctly:

To start the process is usually look at my `CMakeLists.txt` file to figure out what I need.
In this case:
{% highlight bash %}
cd  ParSim/Vaango/src/submodules
ls goolestest
ls json
{% endhighlight %}

The `json` submodule downloads a lot of unnecessary data and will be removed at a future date.

##### Checking needed third party packages #####

Looking at the root-level `CMakeLists.txt` file tells us that the following external packages are
needed to build `Vaango`:

* `boost` : for some MPI and serialization code
* `cmake` : to build the makefiles
* `eigen3` : for some matrix operations
* `gcc` : to build the C++ code
* `gfortran` : to build the Fortran code
* `libxml2` : for XML input/output of `Vaango` format data
* `openmpi` : for MPI code
* `perl` : for Perl scripts
* `zlib` : for compression code

We avoid continuous integration and testing in the build; so `googletests` is not
used even though it is downloaded.

#####  Loading modules #####

A `module avail` command typically lists a large number of potential packages that can
be used.  In our case, we loaded the following modules and environments:
{% highlight bash %}
module load cmake/2.8.10.2
module load gcc/4.9.2
module swap PrgEnv-pgi PrgEnv-gnu
{% endhighlight %}

#####  Installing Boost and Eigen3 #####

I couldn't locate `boost` or `eigen3` on the machine and decided to download and build them:
{% highlight bash %}
wget https://downloads.sourceforge.net/project/boost/boost/1.58.0/boost_1_58_0.tar.bz2
wget http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2
{% endhighlight %}

To build `boost` I did the following:
{% highlight bash %}
tar xvfj boost_1_58_0.tar.bz2
cd boost_1_58_0/
./bootstrap.sh --with-libraries=regex --prefix=<home_directory>/boost-gcc
./b2 install
cd ..
{% endhighlight %}

For `eigen3`, the process was
{% highlight bash %}
tar xvfj 3.3.4.tar.bz2
mkdir eigen3-build
cd eigen3-build
cmake ../eigen-eigen-5a0156e40feb/ -DCMAKE_INSTALL_PREFIX=<home_directory>
{% endhighlight %}

#####  Compiling Vaango #####

I had to modify `CMakeLists.txt` to automatically detect the MPI compilers and settings.
After that, to compile the `Vaango` code, all I needed was

{% highlight bash %}
export CRAYPE_LINK_TYPE=dynamic
mkdir ParSim/Vaango/dbg-gcc
cd ParSim/Vaango/dbg-gcc
cmake ../src -DCMAKE_BUILD_TYPE=Debug \
             -DEIGEN3_INCLUDE_DIR=<home_directory>/include/eigen3 \
             -DBOOST_ROOT=<home_directory>/boost-gcc
make -j4
cd ..
mkdir opt-gcc
cmake ../src -DCMAKE_BUILD_TYPE=Release \
             -DEIGEN3_INCLUDE_DIR=<home_directory>/include/eigen3 \
             -DBOOST_ROOT=<home_directory>/boost-gcc
make -j4
{% endhighlight %}

#### Running the `Vaango` code  on `Copper` ####

To check whether the build produced a working executable, I had to start an interactive 
session with

{% highlight bash %}
qsub -l select=1:ncpus=32:mpiprocs=32 -A <PROJECT_ID> -l walltime=00:30:00 -q debug -X -I
{% endhighlight %}

and then 

{% highlight bash %}
module switch PrgEnv-pgi PrgEnv-gnu
cd ParSim/Vaango
mkdir tests
cd tests
ln -s ../dbg-gcc/StandAlone/vaango vaango_dbg
ln -s ../opt-gcc/StandAlone/vaango vaango_opt
ln -s ../src/StandAlone/inputs/MPM/const_test_hypo.ups .
aprun -n 1 vaango_dbg -mpi ./const_test_hypo.ups
{% endhighlight %}

Larger jobs require the `qsub` queue system and `PBS` scripts.  

#### Building the `Vaango` code on `Excalibur` (Cray XC40) ####

Another Cray system called `Excalibur` is also used to run `Vaango` once in a while.  The
pre-installed packages on this machine vary with time and the following is what had
to be done to get `Vaango` to run on that machine around a year ago.

##### Downloading the code ####

{% highlight bash %}
module load module-git
git clone --recursive https://github.com/bbanerjee/ParSim
{% endhighlight %}

##### Installing cmake ####

We should use at least version 3.2.2 but earlier versions may may with the latest
`Vaango` code. If the build fails, one may need to load some missing modules.
{% highlight bash %}
wget https://cmake.org/files/v3.2/cmake-3.2.2.tar.gz
mkdir localpackages
tar -xvfz cmake-3.2.2.tar.gz
cd cmake-3.2.2
./boostrap --prefix=<home_directory>/localpackages/cmake && make && make install
{% endhighlight %}

##### Installing boost and eigen3 ####

These can be installed in a manner similar to that for Copper.  The
installation directory was chosen to be `localpackages`.

##### Compiling Vaango ####
The cmake script may need the full set of options but typically works with
just the locations of boost and eigen3 provided in the command line, if
the correct environment is chosen  (in our case, gnu).  The process
should be identical to that used in Copper.

The full set of path and library options to cmake when using a local `mpich` installation is given below.

{% highlight bash %}
mkdir opt
cd opt
cmake ../src \
  -DMPI_DIR=/home/banerjee/localpackages/mpich-install \
  -DMPI_C_NO_INTERROGATE:STRING="/home/banerjee/localpackages/mpich-install/bin/mpicc" \
  -DMPI_CXX_NO_INTERROGATE:STRING="/home/banerjee/localpackages/mpich-install/bin/mpicxx" \
  -DMPI_Fortran_NO_INTERROGATE:STRING="/home/banerjee/localpackages/mpich-install/bin/mpifort" \
  -DMPI_C_COMPILER=/home/banerjee/localpackages/mpich-install/bin/mpicc \
  -DMPI_CXX_COMPILER=/home/banerjee/localpackages/mpich-install/bin/mpicxx \
  -DMPI_Fortran_COMPILER=/home/banerjee/localpackages/mpich-install/bin/mpifort \
  -DMPI_C_LIBRARIES:STRING="-lmpi -lmpicxx -L/home/banerjee/localpackages/mpich-install/lib" \
  -DMPI_CXX_LIBRARIES:STRING="-lmpi -lmpicxx -L/home/banerjee/localpackages/mpich-install/lib" \
  -DMPI_Fortran_LIBRARIES:STRING="/home/banerjee/localpackages/mpich-install/lib/libmpifort.so" \
  -DMPI_C_INCLUDE_PATH:STRING="/home/banerjee/localpackages/mpich-install/include"\
  -DMPI_CXX_INCLUDE_PATH:STRING="/home/banerjee/localpackages/mpich-install/include" \
  -DMPI_Fortran_INCLUDE_PATH:STRING="/home/banerjee/localpackages/mpich-install/include" \
  -DMPI_C_LINK_FLAGS:STRING="-L/home/banerjee/localpackages/mpich-install/lib" \
  -DMPI_CXX_LINK_FLAGS:STRING="-L/home/banerjee/localpackages/mpich-install/lib" \
  -DEIGEN3_INCLUDE_DIR=/home/banerjee/localpackages/eigen3-install/include \
  -DBOOST_ROOT=/home/banerjee/localpackages/boost-install
make
{% endhighlight %}

#### Remarks ####

The process of building Vaango on Cray machines has become considerably simpler over time if
implicit codes are not needed.  That is still not true for IBM machines such as BlueGene/Q.



