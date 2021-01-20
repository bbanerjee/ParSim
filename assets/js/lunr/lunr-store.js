var store = [{
        "title": "Blog articles",
        "excerpt":" ","categories": [],
        "tags": [],
        "url": "https://bbanerjee.github.io/docs/articles",
        "teaser": null
      },{
        "title": "Checking the build",
        "excerpt":"To check the build, you should run a few example problems. Create a directory called runs under Vaango mkdir runs Create links to the optimized and debug executables: cd runs ln -s ../dbg/StandAlone/vaango vaango_dbg ln -s ../opt/StandAlone/vaango vaango_dbg Create links to the inputs directory: ln -s ../src/StandAlone/inputs/MPM MPM_inputs Do a...","categories": [],
        "tags": [],
        "url": "https://bbanerjee.github.io/docs/build-check",
        "teaser": null
      },{
        "title": "Building Vaango",
        "excerpt":"Instructions for building Vaango on Ubuntu are given below. Prerequisites Cmake: Compilers: MPI and XML libraries: Other libraries: Building the executables The Vaango repository Out-of-source optimized build Out-of-source debug build Unit tests Clang compiler: Visit build: Compiling the code: Prerequisites Cmake: You will probably need to install Cmake to control...","categories": [],
        "tags": [],
        "url": "https://bbanerjee.github.io/docs/build-instructions",
        "teaser": null
      },{
        "title": "Downloading ParSim",
        "excerpt":"User version Developer version (for Ubuntu) Getting the latest version, adding files etc. User version Create a clone of the ParSim repository git clone https://github.com/bbanerjee/ParSim Check status git status git diff origin/master Just in case the master has changed, update your copy git pull Developer version (for Ubuntu) If you...","categories": [],
        "tags": [],
        "url": "https://bbanerjee.github.io/docs/getting-started",
        "teaser": null
      },{
        "title": "Vaango Manuals",
        "excerpt":"                       The following PDF manuals are based on the most recent version of       Vaango.  Vaango is under continuous development and deployment,       and there are no standardized release cycles.                                       Vaango installation manual                                                      Vaango user manual                                                      Vaango theory manual                                                      Vaango developer manual                                                ","categories": [],
        "tags": [],
        "url": "https://bbanerjee.github.io/docs/manuals",
        "teaser": null
      },{
        "title": "Material and simulation data",
        "excerpt":"Links to some material and simulation data that we have made publicly available. The simulation data consists of input files for various packages. High-rate and high-temperature experiments Our data on steel, aluminum, and copper have been made freely available at https://github.com/bbanerjee/HighRateExptData. Concrete and plastic simulations Input files for the simulation...","categories": [],
        "tags": [],
        "url": "https://bbanerjee.github.io/docs/material-data",
        "teaser": null
      },{
        "title": "Parsim - Particle-based Simulation",
        "excerpt":"Parsim - Particle-based Simulation ParSim is a fork of the Uintah simulation software with a focus on solid mechanics modeling and simulation. Current development is focussed on Vaango, a MPM-CFD code, and GranularSim, a DEM-SPH-Peridynamics code. We are also exploring the development of a user interface for these tools. Documentation...","categories": [],
        "tags": [],
        "url": "https://bbanerjee.github.io/docs/parsim",
        "teaser": null
      },{
        "title": "Tabular plasticity",
        "excerpt":"Shear modulus model Bulk modulus model Yield function Input file formats The prescribed deformation file The bulk-modulus model The yield function The tabular plasticity model is designed to simulate hypoelastic-plasticity using tabulated data for elastic moduli and the yield function for \\(J_2\\) and \\(J_2-I_1\\) models of perfect plasticity. The model...","categories": [],
        "tags": [],
        "url": "https://bbanerjee.github.io/docs/vaango/mpm-materials/tabular-plasticity/",
        "teaser": null
      },{
        "title": "Tabular J2-plasticity with linear elasticity",
        "excerpt":"Coming soon.   ","categories": [],
        "tags": [],
        "url": "https://bbanerjee.github.io/docs/vaango/tutorials/tabular-j2-lin-elas/",
        "teaser": null
      },{
        "title": "Sampling large data sets (Part 1)",
        "excerpt":"The amount of data is large. Before we analyze the whole set of data, we can find trends by examining a smaller set of sensors. A smaller set is also useful if we want to train a model. One of extracting a smaller set from the data is to generate...","categories": ["sampling"],
        "tags": [],
        "url": "https://bbanerjee.github.io/sampling/latin_hypercube_sampling/",
        "teaser": null
      },{
        "title": "Sampling large data sets (Part 2)",
        "excerpt":"Recall that we had found a set of points that samples a cylindrical domain uniformly. However, these points do not match the locations of our sensors and we will have to find the nearest sensor to each sample point. This problem is a version of the Assignment Problem and an...","categories": ["sampling"],
        "tags": [],
        "url": "https://bbanerjee.github.io/sampling/hungarian_algorithm/",
        "teaser": null
      },{
        "title": "Formatting C++ code",
        "excerpt":"Every once in a while I need to recall how clang-format is set up and run on my Ubuntu 16.04 machine. This post is mean to act as a reminder. If I want to change all the files in my repository to a particular format (I prefer the Mozilla style),...","categories": ["C++"],
        "tags": [],
        "url": "https://bbanerjee.github.io/c++/clang-format/",
        "teaser": null
      },{
        "title": "Auto-modernizing C++ code",
        "excerpt":"The Clang static analyzer tools come with a handy interface called clang-tidy. I’ve been try to set up this tool in my cmake toolchain with varying levels of success. In particular I’d like to automate the conversion of old code into a version that uses C++11 and some C++14 constructs....","categories": ["C++"],
        "tags": [],
        "url": "https://bbanerjee.github.io/c++/clang-tidy/",
        "teaser": null
      },{
        "title": "Command pattern for regression testing",
        "excerpt":"A few days ago I had to refactor a 20,000 line class that was being used as a regression tester for discrete element and peridynamics simulations. After some thought, I realized that the easiest way to achieve what I wanted was by using the Command design pattern (minus the undo...","categories": ["C++"],
        "tags": [],
        "url": "https://bbanerjee.github.io/c++/command-pattern/",
        "teaser": null
      },{
        "title": "Input files: Reading XML",
        "excerpt":"Mechanics research codes are typically written by graduate students who aim to get their work done as quickly as possible. These codes are not meant to last beyond the publication of a few related papers. As a result, typical input files have the form 0.000000e+00 -1.000000e+00 0.000000e+00 100 1 50...","categories": ["C++","XML"],
        "tags": [],
        "url": "https://bbanerjee.github.io/c++/xml/xml-input/",
        "teaser": null
      },{
        "title": "Reading JSON in C++",
        "excerpt":"In my previous post on reading XML input files, I discussed how input files can be made more human friendly with XML markup. In some situations, it may be preferable to have JSON format files instead. JSON is particularly useful when the Javascript is used during reading and writing. Though...","categories": ["C++","JSON"],
        "tags": [],
        "url": "https://bbanerjee.github.io/c++/json/json-input/",
        "teaser": null
      },{
        "title": "Writing VTK XML files in C++",
        "excerpt":"I visualize the output of my simulations using either LLNL’s Visit or Kitware’s ParaView. These tools are wonderful for dealing with large datasets and can read a huge variety of file formats. In particular, they are good for remote and visualization. One of the formats that both these tools can...","categories": ["C++","XML","VTK"],
        "tags": [],
        "url": "https://bbanerjee.github.io/c++/xml/vtk/vtk-output/",
        "teaser": null
      },{
        "title": "Particle data in VTK XML",
        "excerpt":"In my previous post on writing VTK output files I described how mesh data can be output in VTK XML format. In this article I will talk about how I output particle data from my simulations. These simulations use a number of techniques depending on the requirements. I use the...","categories": ["C++","XML","VTK"],
        "tags": [],
        "url": "https://bbanerjee.github.io/c++/xml/vtk/vtk-particle-output/",
        "teaser": null
      },{
        "title": "Reading VTK particles in Javascript",
        "excerpt":"In our previous article we discussed the process of calling XML functions in C++ to write data produced by particle simulations in XML format files. These output files can then be viewed in powerful tools such as VisIt or Paraview. Modern browsers have become powerful in recent years and can...","categories": ["Javascript","Typescript","Vue","Vuex","XML"],
        "tags": [],
        "url": "https://bbanerjee.github.io/javascript/typescript/vue/vuex/xml/reading-vtk-particles/",
        "teaser": null
      },{
        "title": "Saving particle data in a Vuex store",
        "excerpt":"Introduction Typescript typings for Vue and Vuex The store The particle module The particle state The particle getters The particle mutations The particle actions Remarks Introduction In the first part of this series on scientific visualization in Javascript, we saw how data in VTK unstructured grid format produced by C++...","categories": ["Javascript","Typescript","Vue","Vuex"],
        "tags": [],
        "url": "https://bbanerjee.github.io/javascript/typescript/vue/vuex/vuex-store/",
        "teaser": null
      },{
        "title": "Plotting VTK particles with Three.js",
        "excerpt":"Introduction A Vuex store for graphics Interacting with Vue Conclusion Introduction In Part 2 of this series, we showed how the particle data are saved into a Vuex store. Now we are ready to visualize the data. In this article, we will discuss how three.js can be used to display...","categories": ["Javascript","Typescript","Threejs","Vue"],
        "tags": [],
        "url": "https://bbanerjee.github.io/javascript/typescript/threejs/vue/vtk-threejs/",
        "teaser": null
      },{
        "title": "Setting up the Three.js renderer",
        "excerpt":"Introduction The graphics panel The renderer template The renderer implementation Remarks Introduction The previous articles in this series were about: Part 1: Reading VTK format particles with Javascript in a browser Part 2: Saving the read-in particle data in a Vuex store Part 3: Initialization of a store and the...","categories": ["Javascript","Typescript","Threejs","Vue"],
        "tags": [],
        "url": "https://bbanerjee.github.io/javascript/typescript/threejs/vue/vtk-threejs-renderer/",
        "teaser": null
      },{
        "title": "Setting up the Three.js scene and camera",
        "excerpt":"Introduction Flash-back: The graphics panel The three-scene component The three-camera component Remarks Introduction In the previous articles in this series we talked about: Part 1: Reading VTK format particles with Javascript in a browser Part 2: Saving the read-in particle data in a Vuex store Part 3: Initialization of a...","categories": ["Javascript","Typescript","Threejs","Vue"],
        "tags": [],
        "url": "https://bbanerjee.github.io/javascript/typescript/threejs/vue/vtk-threejs-camera/",
        "teaser": null
      },{
        "title": "Setting up the Three.js ellipsoids",
        "excerpt":"Introduction The ellipsoid particles component Creating the ellipsoid particles Remarks Introduction In the previous articles in this series we talked about: Part 1: Reading VTK format particles with Javascript in a browser Part 2: Saving the read-in particle data in a Vuex store Part 3: Initialization of a store and...","categories": ["Javascript","Typescript","Threejs","Vue"],
        "tags": [],
        "url": "https://bbanerjee.github.io/javascript/typescript/threejs/vue/vtk-threejs-ellipsoid/",
        "teaser": null
      },{
        "title": "Plotting particles with vtk.js",
        "excerpt":"Introduction Registering VTK components Creating Vuex stores for VTK data The VTK graphics panel The VTK renderer The VTK particles component Creating the particles Remarks Introduction In the previous articles in this series we talked about: Part 1: Reading VTK format particles with Javascript in a browser Part 2: Saving...","categories": ["Javascript","Typescript","vtkjs","Threejs","Vue"],
        "tags": [],
        "url": "https://bbanerjee.github.io/javascript/typescript/vtkjs/threejs/vue/vtkjs-ellipsoid/",
        "teaser": null
      },{
        "title": "XML format for particle input files",
        "excerpt":"Introduction R script converter from CSV to XML The output ASCII XML file The output compressed base64 XML file Remarks Introduction Let us now change tack slightly and return to an issue I had talked about earlier in Input files: reading XML. Typical input files in research simulation codes cannot...","categories": ["R","XML"],
        "tags": [],
        "url": "https://bbanerjee.github.io/r/xml/xml-particle-input-file/",
        "teaser": null
      },{
        "title": "Plane stress Drucker-Prager return algorithm",
        "excerpt":"Introduction Plane stress elasticity Drucker-Prager plasticity Plane stress Drucker-Prager yield function Plane stress Drucker-Prager flow rule Remarks Introduction Recently I’ve encountered questions on how the radial return algorithm works when applied to plane stress plasticity. There seems to be some confusion about the application of the plane strain constraint. Here’s...","categories": ["Mechanics","Plasticity","Algorithm"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mechanics/plasticity/algorithm/plane-stress-drucker-prager/",
        "teaser": null
      },{
        "title": "The plane stress return algorithm",
        "excerpt":"Introduction Review of 3D plasticity Plane stress plasticity Forward Euler Return algorithm Remarks Introduction In the previous part of this discussion, I derived plane stress expressions for linear elasticity, the Drucker-Prager yield function, and the associated flow rule. Let us now review the approach used for finding the parameter \\(\\dot{\\lambda}\\)...","categories": ["Mechanics","Plasticity","Algorithm"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mechanics/plasticity/algorithm/plane-stress-drucker-prager-return/",
        "teaser": null
      },{
        "title": "Plane stress return: Spectral decomposition",
        "excerpt":"Introduction Spectral decomposition of \\(\\mathbf{C}\\) Spectral decomposition of \\(\\mathbf{P}\\) Product of \\(\\mathbf{P}\\) and \\(\\mathbf{C}\\) Remarks Introduction The previous parts of this series dealt with: Part 1: Background of 3D and plane stress Drucker-Prager plasticity Part 2: A forward Euler return algorithm for plane stress plasticity Let us now take a...","categories": ["Mechanics","Plasticity","Algorithm"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mechanics/plasticity/algorithm/plane-stress-drucker-prager-return-part-3/",
        "teaser": null
      },{
        "title": "Plane stress forward Euler Drucker-Prager",
        "excerpt":"Introduction Finding \\(\\Delta\\lambda\\) An attempt at simplification Remarks Introduction In Part 2 of this series, we saw that a forward Euler return algorithm leads to the discretized equations $$ \\boldsymbol{\\sigma}_{n+1} = \\boldsymbol{\\sigma}_{n+1}^{\\text{trial}} - \\Delta\\lambda \\,\\mathbf{C}\\, \\boldsymbol{n}_{n} $$ where $$ \\boldsymbol{n}_n = \\begin{bmatrix} n^n_{11} \\\\ n^n_{22} \\\\ n^n_{12} \\end{bmatrix} = \\left(...","categories": ["Mechanics","Plasticity","Algorithm"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mechanics/plasticity/algorithm/plane-stress-drucker-prager-return-part-4/",
        "teaser": null
      },{
        "title": "Forward vs. Backward Euler: Plane stress plasticity",
        "excerpt":"Introduction Forward difference for stress rate Forward and Backward Euler for stress rate Forward and Backward Euler for flow rule Forward/Backward Euler stress updates Does the choice of Forward/Backward Euler matter? Accuracy Stability Optimization Remarks Introduction One of the main points of divergence of many implementations of plastic return algorithms...","categories": ["Mechanics","Plasticity","Algorithm"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mechanics/plasticity/algorithm/plane-stress-drucker-prager-return-part-5/",
        "teaser": null
      },{
        "title": "Nonlinear programming and closest point return plasticity",
        "excerpt":"Introduction Background Primal form The Lagrangian Dual function Dual form Karush-Kuhn-Tucker optimality conditions Similarity with plasticity Closest point return Remarks Introduction In Part 5 I briefly hinted at the closest-point return algorithm. The ideas behind this were made rigorous in the mid-to-late 1980s by a group of researchers influenced by...","categories": ["Mechanics","Plasticity","Algorithm"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mechanics/plasticity/algorithm/primal-dual-closest-point-return/",
        "teaser": null
      },{
        "title": "Exploring closest point return plasticity",
        "excerpt":"Introduction Eigendecompositions in linear elasticity The transformed space for isotropic linear elasticity The Lode invariants and the Lode basis The transformed stress tensor Remarks Introduction In Part 6, I explained why a backward Euler stress update and a closest point return from the trial stress to the yield surface are...","categories": ["Mechanics","Plasticity","Algorithm"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mechanics/plasticity/algorithm/closest-point-return/",
        "teaser": null
      },{
        "title": "Geometric closest point return algorithm",
        "excerpt":"Introduction The Arena yield function The non-hardening return algorithm The closest point algorithm An animation of the closest point algorithm Remarks Introduction In Part 7, we saw that for isotropic elastic materials and perfect associated plasticity, the trial stress and the actual stress are at the shortest distance from each...","categories": ["Mechanics","Plasticity","Algorithm"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mechanics/plasticity/algorithm/geometric-closest-point-return/",
        "teaser": null
      },{
        "title": "Creating an animation with d3.js",
        "excerpt":"Introduction Input data The HTML file The Javascript code Generating points on the yield surface Drawing the yield surface and closest-point projection Creating the SVG Creating the canvas group Creating map from real to canvas coordinates Creating SVG polyline generator Creating the axes and the yield surface Creating groups for...","categories": ["Javascript","D3JS"],
        "tags": [],
        "url": "https://bbanerjee.github.io/javascript/d3js/d3-animation-closest-point-return/",
        "teaser": null
      },{
        "title": "Reading XML files containing gzipped data in C++",
        "excerpt":"Introduction Recap The header file The implementation The read function The readParticleValues templated function The decodeAndUncompress templated function The convert&lt;T&gt; template specializations Remarks Introduction We saw how to create an XML file containing compressed particle data in the article “XML format for particle input files”. Let us now explore how...","categories": ["C++","XML","gzip"],
        "tags": [],
        "url": "https://bbanerjee.github.io/c++/xml/gzip/reading-xml-with-gzipped-data/",
        "teaser": null
      },{
        "title": "The CFL condition for explicit discrete element methods:1",
        "excerpt":"Introduction The CFL condition Remarks Introduction Solutions of hyperbolic partial differential equations using explicit numerical methods need a means of limiting the timestep so that the solution is stable. A criterion that is usually used to constrain the step size is the Courant–Friedrichs–Lewy condition. Let us first explore what the...","categories": ["DEM"],
        "tags": [],
        "url": "https://bbanerjee.github.io/dem/CFL-condition-discrete-elements/",
        "teaser": null
      },{
        "title": "The CFL condition for explicit discrete element methods:2",
        "excerpt":"Introduction The von Neumann approach to stability analysis Remarks Introduction In the first part of this article, we looked at the one-dimensional linear second-order wave equation $$ \\frac{\\partial^2 u}{\\partial t^2} - \\frac{\\partial^2 u}{\\partial x^2} = 0 $$ We saw that for a discretization with grid size \\(\\Delta x\\) and timestep...","categories": ["DEM"],
        "tags": [],
        "url": "https://bbanerjee.github.io/dem/CFL-condition-discrete-elements-part-2/",
        "teaser": null
      },{
        "title": "The CFL condition for explicit discrete element methods:3",
        "excerpt":"Introduction Euler’s laws of motion Discretization of Euler equations Stability of central difference scheme One-dimensional unforced system Stability using the Hurwitz matrix approach Remarks Introduction In Part 1 of this article, we revisited the CFL condition and in Part 2 we showed how the CFL condition and the von Neumann...","categories": ["DEM"],
        "tags": [],
        "url": "https://bbanerjee.github.io/dem/CFL-condition-discrete-elements-part-3/",
        "teaser": null
      },{
        "title": "The CFL condition for explicit discrete element methods:4",
        "excerpt":"Introduction One-dimensional equations for impact Stability of central difference for one-dimensional impact Remarks Introduction In Part 3 of this article I discussed the approach where the equations for a system of rigid bodies are approximated by a spring-mass system. The numerical stability conditions of that system are then taken to...","categories": ["DEM"],
        "tags": [],
        "url": "https://bbanerjee.github.io/dem/CFL-condition-discrete-elements-part-4/",
        "teaser": null
      },{
        "title": "Are stresses tensile or compressive during rigid body rotation?",
        "excerpt":"Introduction The question The finite element solution Remarks Introduction Recently, I came across a question in StackExchange that pointed out that some books on continuum mechanics suggest that an element will increase in size when rotated if a small strain approximation is used in a finite element simulation. This issue...","categories": ["FEM"],
        "tags": [],
        "url": "https://bbanerjee.github.io/fem/rigid-body-rotation-small-strain/",
        "teaser": null
      },{
        "title": "Material and spatial incremental constitutive equations",
        "excerpt":"The question Instantaneous moduli for PK-2 stress and Green strain Instantaneous moduli for Kirchhoff stress The reason for the inconsistency The question A colleague asked a question on objectivity a few days ago that had me going back to Ray Ogden’s book on nonlinear elastic deformations. The question was on...","categories": ["Mechanics"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mechanics/incremental-objectivity/",
        "teaser": null
      },{
        "title": "The ARENA model for partially saturated soils",
        "excerpt":"Introduction Some predictions from ARENA Remarks Introduction I developed the ARENA model for partially saturated soils last year (2016). Before that we had tried using a typical Drucker-Prager with cap model and then a modified Cam-Clay model but could not get these models to either represent observed experimental data or...","categories": ["Mechanics"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mechanics/Arena-model-paper/",
        "teaser": null
      },{
        "title": "The difference between the spin and angular velocity tensors",
        "excerpt":"Introduction The spin tensor The Green-Naghdi objective rate Relation between spin and angular velocity Remarks Introduction Several years ago I took created a Wikpedia article on objective stress rates based on lecture notes in nonlinear finite elements for a class that I taught at the University of Utah. In this...","categories": ["Mechanics"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mechanics/spin-tensor-and-Green-Naghdi-rate/",
        "teaser": null
      },{
        "title": "Can the Larsen-C ice shelf failure be predicted with Peridynamics?",
        "excerpt":"Introduction Issues to be resolved before a simulation Remarks Introduction On July 7, 2017, the Project Midas group released a couple of plots of material velocities and interferograms showing the evolution of a large rift the Larsen-C ice shelf in Antarctica. The crack grew 11 miles in a few days...","categories": ["Fracture"],
        "tags": [],
        "url": "https://bbanerjee.github.io/fracture/Larsen-B-and-peridynamics/",
        "teaser": null
      },{
        "title": "Parallel domain decomposition for particle methods: Part 1",
        "excerpt":"Introduction Creating and scattering particles MPI implementation MPI setup The scatter operation Remarks Introduction For parallel particle codes that have to be written quickly (while retaining flexibility), the task-based parallelism approach doesn’t always work well. The usual approach that is taken in those situations is some sort of domain decomposition...","categories": ["MPI","C++"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mpi/c++/parallel-domain-decomposition-part-1/",
        "teaser": null
      },{
        "title": "Parallel domain decomposition for particle methods: Part 2",
        "excerpt":"Introduction Exchanging particles between processes MPI implementation PatchNeighborComm struct Patch struct The particle exchange function Remarks Introduction The previous article in this series discussed the scatter operation for moving particles to various processes. In this second part of the series we will discuss a commonly used method of communicating information...","categories": ["MPI","C++"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mpi/c++/parallel-domain-decomposition-part-2/",
        "teaser": null
      },{
        "title": "Parallel domain decomposition for particle methods: Part 3",
        "excerpt":"Introduction Plimpton’s scheme for exchanging particles MPI implementation Patch struct The particle exchange function Remarks Introduction In Part 2 of this series we showed the direct way of communicating ghost particles between patches. That approach requires 26 communication steps per patch in three-dimensions. In this article we discuss the approach...","categories": ["MPI","C++"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mpi/c++/parallel-domain-decomposition-part-3/",
        "teaser": null
      },{
        "title": "Parallel domain decomposition for particle methods: Part 4",
        "excerpt":"Introduction Plimpton’s scheme for migrating particles MPI implementation Patch struct Remarks Introduction The Plimpton scheme of communicating ghost information between patches was described in Part 3 of this series. Let us now see how a similar approach can be used to migrate particles that have moved across patches. In the...","categories": ["MPI","C++"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mpi/c++/parallel-domain-decomposition-part-4/",
        "teaser": null
      },{
        "title": "Unit testing with MPI, googletest, and cmake",
        "excerpt":"Introduction Installing googletest Making sure cmake finds and compiles googletest Adding local unit tests The CMakeLists.txt file in UnitTests The actual test C++ code The MPI test environment class The main test function The actual test Caveat The output from make Remarks Introduction In this article we take a short...","categories": ["MPI","C++"],
        "tags": [],
        "url": "https://bbanerjee.github.io/mpi/c++/mpi-unit-testing-googletests-cmake/",
        "teaser": null
      },{
        "title": "Compiling and running the MPM code Vaango on a Cray",
        "excerpt":"Introduction Authentication Building the Vaango code on Copper (Cray XE6m) Downloading the code Checking needed third party packages Loading modules Installing Boost and Eigen3 Compiling Vaango Running the Vaango code on Copper Building the Vaango code on Excalibur (Cray XC40) Downloading the code Installing cmake Installing boost and eigen3 Compiling...","categories": ["Vaango"],
        "tags": [],
        "url": "https://bbanerjee.github.io/vaango/compiling-vaango-on-a-cray/",
        "teaser": null
      },]
