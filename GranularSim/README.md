DEM + PD + SPH

# License #
This code has been derived from Tahoe.  The license terms for Tahoe can be found in
UCB.license.  The modified code has an addition license that can be found in Parresia.license.

# Caveats/Special Requirements #
This code uses cppcodec to base64 conversions.  The unit tests in this code use googletests.
Both of these dependencies are included as submodules in this git repository.

To make sure that these dependencies are checked out when cloning this code, use

git clone --recursive



How to use paraEllip3d

1. compile

   1) setup compile environment, e.g.:
      module load openmpi-1.6.4-gcc-4.6.4

   2) setup boost in makefile, e.g.:
      BOOST_ROOT = /usr/local/boost-1.53.0-openmpi-1.6.4-gcc-4.6.4
      Note: run environment must conform, for example, 
            Torque/PBS qsub script should have the line:
            module load openmpi-1.6.4-gcc-4.6.4

   3) make

   4) copy paraEllip3d to your simulation directory

2. run

   1) 1-process-1-thread mode

      ./paraEllip3d input.txt or
      mpirun -np 1 ./paraEllip3d input.txt

      In input.txt:
       mpiProcX  1
       mpiProcY  1
       mpiProcZ  1
       ompThreads 1

   2) 1-process-multi-thread OpenMP mode

      ./paraEllip3d input.txt or
      mpirun -np 1 ./paraEllip3d input.txt

      In input.txt:
       mpiProcX  1
       mpiProcY  1
       mpiProcZ  1
       ompThreads 12

   3) MPI-only mode (recommended on a single workstation)

      mpirun -np 64 ./paraEllip3d input.txt

      In input.txt:
       mpiProcX  4
       mpiProcY  4
       mpiProcZ  4
       ompThreads 1

      It must satisfy np = mpiProcX * mpiProcY * mpiProcZ,
      i.e., 64 = 4 * 4 * 4, where np obviously means number 
      of processes. And ompThreads = 1 must be specified.

   4) MPI/OpenMP hybrid mode (recommended on clusters)

      mpirun -np 64 -npernode 1 ./paraEllip3d input.txt

      In input.txt:
       mpiProcX  4
       mpiProcY  4
       mpiProcZ  4
       ompThreads 12

      It must satisfy np = mpiProcX * mpiProcY * mpiProcZ,
      i.e., 64 = 4 * 4 * 4, where np not only means number of 
      processes, but also means number of nodes.

      "-npernode 1" means each node runs only 1 process, and 
      that process contains multiple threads specified by 
      ompThreads = n (n > 1).

3. job script examples for Torque/PBS

   1) On soilblast

    !/bin/sh

    #PBS -m abe
    #PBS -M your_email_address
    #PBS -j oe

    #PBS -N benchmark
    #PBS -l nodes=1:ppn=12
    #PBS -l walltime=24:00:00

    module purge
    module load openmpi-1.6.4-gcc-4.6.4

    cd $PBS_O_WORKDIR
    mpirun -np 12 ./paraEllip3d input.txt
    
   2) On Janus(it is a good practice to compile and run 
      in the same environment, i.e., use job script for
      compile.)

    #!/bin/sh

    #PBS -m abe
    #PBS -M your_email_address
    #PBS -j oe

    #PBS -N benchmark
    #PBS -l nodes=8:ppn=12
    #PBS -l walltime=24:00:00
    #PBS -q janus-small
    #PBS -A CEAE00000001

    . /curc/tools/utils/dkinit
    reuse Torque
    reuse Moab
    reuse .openmpi-1.4.3_gcc-4.5.2_torque-2.5.8_ib

    cd $PBS_O_WORKDIR
    mpirun -np 8 -npernode 1 ./paraEllip3d input.txt


#-----------------------------------------------------------------------
# Formatting the code using clang-format
#-----------------------------------------------------------------------
Run the following:

    find . -name "*.cpp" -exec clang-format -i -style Mozilla {} \;
    find . -name "*.h" -exec clang-format -i -style Mozilla {} \;


