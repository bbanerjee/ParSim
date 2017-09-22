---
layout: post-right-sidebar
title:  "Unit testing with MPI, googletest, and cmake"
subheadline: "Biswajit Banerjee"
description: "Continuous integration with make"
date:  2017-08-11 10:30:00
categories:
    - MPI
    - C++
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

#### Introduction ####
In this article we take a short detour into the problem of continuous unit testing
of code that contains MPI calls and use either [mpich](https://www.mpich.org/)
or [Open MPI](https://www.open-mpi.org/).  I have recently moved from
[CTest-based testing](https://cmake.org/Wiki/CMake/Testing_With_CTest) to
a combination of [CMake](https://cmake.org/) and 
[googletest](https://github.com/google/googletest).  The reason for the shift
is the convenience `googletest`
provides.  However, the `googletest`
[primer](https://github.com/google/googletest/blob/master/googletest/docs/Primer.md)
and [advanced manual](https://github.com/google/googletest/blob/master/googletest/docs/AdvancedGuide.md)
do not contain many examples and can be cryptic at times.  I hope this
article will provide pointers to those who run into a roadblock when testing MPI
applications with `googletest`.

#### Installing googletest ####
I typically install `googletest` as a [submodule](https://git-scm.com/docs/git-submodule)
in my git repository.  For example,

{% highlight bash %}
git submodule add git@github.com:google/googletest.git ./googletest
{% endhighlight %}

There are several [caveats](https://chrisjean.com/git-submodules-adding-using-removing-and-updating/) 
when using git submodules.  For our purposes, we only have to remember that when we clone 
the repository we have to use

{% highlight bash %}
git clone --recursive
{% endhighlight %}

Other tips can be found in a nicely condensed form in [Sentheon's blog](https://sentheon.com/blog/git-cheat-sheet.html#working-with-submodules).

#### Making sure cmake finds and compiles googletest ####
To make sure that `googletest` can be found and is built during the compile process,
I add the following to by root `CMakeLists.txt` file:

{% highlight cmake %}
if (USE_CLANG)
  set(CMAKE_CXX_COMPILER "/usr/local/bin/clang++")
endif ()
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
add_subdirectory(googletest)
{% endhighlight %}

The first line is needed because I have not been successful in passing the `clang` compiler
name to `googletest` without specifying the full path.  I'm sure one can use a more general
approach, but I haven't felt the need to spend the time trying to figure out a better way.
The `add_subdirectory`	command is all that is needed for `cmake` to compile `googletest` and
produce static libraries.

#### Adding local unit tests ####
Next I add a `UnitTests` directory in my directory of interest and modify the local
`CMakeLists.txt` to be able to find the `UnitTests` directory.  For example,

{% highlight cmake %}
SET(SPH_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/SmoothParticleHydro.cpp
)
SET(ELLIP3D_SRCS
  ${ELLIP3D_SRCS}
  ${SPH_SRCS}
  PARENT_SCOPE
)
add_subdirectory(UnitTests)
{% endhighlight %}

In this case I am going to test my
[Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) code.

#### The CMakeLists.txt file in UnitTests ####
Now we are finally really to add our unit tests to the build chain.  The `CMakeLists.txt`
file in the `UnitTests` directory has two sections (which can be simplified if you
are so inclined).

In the first section, we find the location of the `googletest` headers and libraries:

{% highlight cmake %}
set(GTEST_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/googletest/googletest/include")
set(GTEST_LIB gtest_main gtest)
include_directories(${GTEST_INCLUDE_DIR})
{% endhighlight %}

Note that the two `googletest` libraries that we use are `gtest_main` and `gtest`.

In the second section, we add the actual test that requires `MPI`.  Once again, these details
can be abstracted into a `cmake` function if you so desire.

{% highlight cmake %}
add_executable(testSPHParticleScatter testSPHParticleScatter.cpp)
target_link_libraries(testSPHParticleScatter
  ${GTEST_LIB}
  ellip3D_lib
  ${MPI_LIBRARY}
  ....
)
set(UNIT_TEST testSPHParticleScatter)
set(MPI_COMMAND mpirun -np 2 ${UNIT_TEST})
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${MPI_COMMAND}
{% endhighlight %}

The `add_executable` line identifies the unit test program `testSPHParticleScatter.cpp`
which tests the scatter operation between two MPI processes.

The `target_link_libraries` lists the libraries that are needed: the `googletest`
libraries (`GTEST_LIB`) and our coupled DEM-SPH code library (`ellip3D_lib`) 

Then we set up a command to run (`MPI_COMMAND`) and make it use `mpirun`.  You
can generalize this if you want.

Finally we, add the `add_custom_command` line that tells `cmake` to run the `MPI_COMMAND`
after the build is complete (`POST_BUILD`).

#### The actual test C++ code ####
Now that the build system has been configured, we just write our unit test
`testSPHParticleScatter.cpp`.

First, we include the required headers:

{% highlight cpp %}
#include <SmoothParticleHydro/SmoothParticleHydro.h>
#include <gtest/gtest.h>
#include <boost/mpi.hpp>
{% endhighlight %}

Note that we are using the [Boost MPI](http://www.boost.org/doc/libs/1_64_0/doc/html/mpi.html) wrappers.

##### The MPI test environment class #####
But we cannot use the `boost::mpi::environment` call to set up `MPI` (because the `environment`
object is deleted before `googletest`s are run).  Instead, we have to set up a custom environment
to run our `MPI` tests by creating a `MPIEnvironment` class that extends the `::testing::Environment`
class provided by `googletest`.

{% highlight cpp %}
class MPIEnvironment : public ::testing::Environment
{
public:
  virtual void SetUp() {
    char** argv;
    int argc = 0;
    int mpiError = MPI_Init(&argc, &argv);
    ASSERT_FALSE(mpiError);
  }
  virtual void TearDown() {
    int mpiError = MPI_Finalize();
    ASSERT_FALSE(mpiError);
  }
  virtual ~MPIEnvironment() {}
};

{% endhighlight %}

Here, the `SetUp` function calls `MPI_Init` and sets up the environment while the `TearDown` function
calls 'MPI_Finalize`.  All tests are performed when the environment is active.

##### The `main` test function #####
The `main` function in typically not needed in standard `googletest` tests and is generate
by some internal magic.  However, we do need a `main` function in `testSPHParticleScatter.cpp`
because we are using `mpirun` to run the test.

{% highlight cpp %}
int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
  return RUN_ALL_TESTS();
}
{% endhighlight %}

The MPI specific environment is created by the `AddGlobalTestEnvironment` function to which we
pass a `MPIEnvironment` object.  The object is deleted internally by `googletest`.

The tests in `testSPHParticleScatter.cpp` are then run using `RUN_ALL_TESTS()`.  Note that this
function returns from `main` and, therefore, and non-googletest environment objects created in
main (e.g., with boost::mpi::environment) are deleted before the tests are run.

##### The actual test #####
Finally, we add an actual test to `testSPHParticleScatter.cpp` as follows:

{% highlight cpp %}
TEST(SPHParticleScatterTest, scatter)
{
  // Set up communicator
  boost::mpi::communicator boostWorld;
  // Set up SPH object
  SmoothParticleHydro sph;
  sph.setCommunicator(boostWorld);
  // ... Some code to create particles
  // ....
  sph.setParticles(particles);
  sph.scatterSPHParticle(domain, ghostWidth, domainBuffer);
  if (boostWorld.rank() == 0) {
    EXPECT_EQ(sph.getSPHParticleVec().size(), 100);
  } else {
    EXPECT_EQ(sph.getSPHParticleVec().size(), 200);
  }
}
{% endhighlight %}

The test can be written in the standard way and numerous examples can be found on the web.
Of course, we have to be careful about keeping in mind that the test will be run
on two processes in this particular case.

##### Caveat #####
The `Boost` MPI environment created by
[`boost::mpi::environment`](http://www.boost.org/doc/libs/1_50_0/boost/mpi/environment.hpp)
allows MPI calls to
fail without throwing an exception.  For example, in my `Patch` code discussed in
an earlier article, I have

{% highlight cpp %}
MPI_Cart_rank(cartComm, neighborCoords.data(), &neighborRank);
{% endhighlight %}

I deliberately send invalid `neighborCoords` to this function to find out if a patch
is a boundary patch.  This does not create a problem when I use the `boost::mpi::environment`
set up.  But when setting up the environment explicitly for the unit tests, I have
to add

{% highlight cpp %}
MPI_Errhandler_set(cartComm, MPI_ERRORS_RETURN);
{% endhighlight %}

before I call `MPI_Cart_rank` to make sure I don't get errors when I use invalid
values of `neighborCoords`.

#### Remarks ####
I hope this article has been of use to you.  Our series on communication between patches
will continue when I get some free time.


If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>


