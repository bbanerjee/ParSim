Verification using the Aldridge solution
========================================
There are three steps in the verification process:

1.  [Generating the Aldridge solution](#generating-the-aldridge-solution).
2.  [Generating the equivalent MPM solution](#generating-the-mpm-solution).
3.  [Comparing the Aldridge solution with the MPM solution](#comparing-the-aldridge-and-mpm-solutions).

Generating the Aldridge solution
--------------------------------
The suggested procedure is as follows:

### Compile the code and set up directories for code execution
1.  Compile the `explode_for.f` file in the `AldridgeSphericalSourceVerificationSolver` directory using 
    the `makefile` in that directory.  This is an in-source build.
2.  Create a directory in `ParSim/Vaango/runs` to run the verification simulations.  We will use the name
    `Verification` as an example.
3.  Go to `Parsim/Vaango/runs/Verification` and create a directory `AldridgeRuns`.
4.  Create a link to `explode_for.exe` in `Parsim/Vaango/runs/Verification/AldridgeRuns`.
```sh
  ln -s src/<path>/AldridgeSphericalSourceVerificationSolver/explode_for.exe explode_for.exe
```
5.  Copy `runAldridge.sh` to `Parsim/Vaango/runs/Verification/AldridgeRuns`. 
6.  Copy `writeAldridgeInput.R` and `plotAndWriteLoadCurve.R` to `Parsim/Vaango/runs/Verification`.

### Create input files 
1.  Go to `Parsim/Vaango/runs/Verification`
2.  Run the R script `writeAldridgeInput.R` using 
```sh
  Rscript writeAldridgeInput.R
```
This creates two files: `input_param_file_hires.in` and `input_param_file_lores.in` which correspond to
the slightly different locations of MPM material points and grid points for two different grid resolutions.
3.  In addition to the input parameter files, we also need the applied load curve in a format that 
can be understood by `explode_for.exe`.  We will use the Vaango UPS file `CentrifugeLoadCurveLow.xml` in
this example.
4. Run the R script `plotAndWriteLoadCurve.R`.  It is preferred that you run this from a R command prompt
using
```R
   source("plotAndWriteLoadCurve.R")
```
because some plots are produced during the process.  The output file generated by this code will be
`CentrifugeLoadCurveLow.xml.resampled`.

### Run the Aldridge executable
1. Go to `Parsim/Vaango/runs/Verification/AldridgeRuns`.
2. Run the shell script
```sh
  sh runAldridge.sh
```
3. It should take a fraction of a second to runs the two cases.  The output files are 
   `output_data_hires.out` and `output_data_lores.out`.  These contain the Aldridge solutions for the
   two cases.

Generating the MPM solution
--------------------------------
1. Go back to `Parsim/Vaango/runs/Verification`.
2. The input files for the verification tests in this directory are
 * Spherical_source_verification_gimp.ups
 * Spherical_source_verification_gimp_quarter.ups
 * Spherical_source_verification_gimp_quarter_hires.ups
 * Spherical_source_verification_gimp_quarter_hires_cpdi.ups
 * Spherical_source_verification_gimp_quarter_hires_cpdi_nocbdi.ups
 * Spherical_source_verification_gimp_quarter_hires_cpdi_nocbdi_nodamp.ups
3. These input file depend on
 * Arena_PhysicalProperties.xml
 * Arena_Arenisca3Parameters.xml
 * CentrifugeLoadCurveLow.xml
4. To run these tests use a command of the form
```sh
  nohup mpirun -np 12 ParSim/Vaango/opt/StandAlone/vaango Spherical_source_verification_gimp.ups > Spherical_source_verification_gimp.out 2>& 1 &
```
5. After running for some time, an output file `Spherical_source_verification_gimp.uda.000` will be 
   produced.
6. Run the shell script `extract_fullspace.sh` to extract the relevant data from these MPM runs
   after uncommenting the appropriate lines.  The associated Python scripts are
   * extract_fullspace_acc.py
   * extract_fullspace_partdata.py
   * extract_fullspace_press.py
7. These runs will produce text files in a `DataAndPlots` subdirectory.

Comparing the Aldridge and MPM solutions
----------------------------------------
1.  Copy the relevant output files from the `DataAndPlots` directory to `Parsim/Vaango/runs/Verification`.
    You can find the names of these files by examining the R script `plotAldridge.R` in
    the directory `Parsim/Vaango/runs/Verification/AldridgeRuns`.
2.  Go to `Parsim/Vaango/runs/Verification/AldridgeRuns`.
3.  Run `plotAldridge.R` after starting R.  Use
```R
  source("plotAldridge.R")
```
4.  The run will produce comparisons between the Aldridge solution and the MPM solution in the form
    of PDF plots.




