Verification using the 1-D wave equation
========================================
There are three steps in the verification process:

1.  [Generating the 1D solution](#generating-the-1d-solution).
2.  [Generating the equivalent MPM solution](#generating-the-mpm-solution).
3.  [Comparing the 1D solution with the MPM solution](#comparing-the-1d-and-mpm-solutions).

Generating the 1D solution
--------------------------------
Solutions to the 1D wave equation for elastic bodies can be found in 
http://imechanica.org/files/awave.pdf

The most convenient solution to use is the one where a velocity vs. time curve is
applied at one end of the bar.  

Generating the MPM solution
--------------------------------
The procedure for generating the MPM solution is:
1. Copy one of the input files listed below to `Parsim/Vaango/runs` or in a directory 
   of your choice that you want your simulation results to be saved in.
2. Create links to the `vaango` and `lineextract` executables in that directory.
3. The input files for the 1D verification tests are listed below.  These can be used
   to compare the effect of different simulation options on the results.
  * Option 1: Velocity BC applied at the domain boundary
    * OneD_velocityBC_gimp.ups
    * OneD_velocityBC_gimp_momform.ups
    * OneD_velocityBC_gimp_damped.ups
    * OneD_velocityBC_gimp_damped_uintah.ups
  * Option 2: Traction BC applied at one end of the bar.  These use `OneD_tractionBC.xml`.
    * OneD_tractionBC_gimp.ups
    * OneD_tractionBC_gimp_damped.ups
  * Option 3: Velocity BC applied through impact
    * OneD_impact_velocityBC_gimp.ups
      * Needs the file `OneD_impact_velocityBC.txt` for specifying the velocity
    * OneD_impact_velocityBC_gimp_damped.ups
      * Needs the same file as above for specifying velocity
    * OneD_impact_velocityBC_hat_gimp_damped.ups
      * Needs the file `OneD_impact_velocityBC_hat.txt` for specifying the velocity
  * Option 4: Pressure BC derived from impact reaction forces and applied at one end.  These use
    `RigidReactionForceHat.xml`.
    * OneD_pressureBC_from_velBC_gimp_damped.ups
    * OneD_pressureBC_from_velBC_hat_gimp_damped_arenisca_elastic_lores.ups
    * OneD_pressureBC_from_velBC_hat_gimp_damped_arenisca_elastic.ups
    * OneD_pressureBC_from_velBC_hat_gimp_damped.ups
    * OneD_pressureBC_from_velBC_hat_gimp_undamped_arenisca_elastic_lores.ups
    * OneD_pressureBC_from_velBC_hat_gimp_undamped_arenisca_elastic.ups
  * Option 5: Pressure BC from centrifuge experimental data
    * OneD_pressureBC_explosion_cpdi_damped_lores.ups
    * OneD_pressureBC_explosion_gimp_damped_hires.ups
    * OneD_pressureBC_explosion_gimp_damped_lores_cfl.ups
    * OneD_pressureBC_explosion_gimp_damped_lores.ups
    * OneD_pressureBC_explosion_gimp_damped.ups
    * OneD_pressureBC_explosion_gimp_undamped_lores.ups
    * OneD_pressureBC_explosion_gimp_undamped.ups
4. These input files may also need 
  * Arena_PhysicalProperties.xml
  * ParamAreniscaElastic.xml
  * CentrifugeLoadCurveLow.xml
5. To run these tests use a command of the form
```sh
  vaango OneD_pressureBC_from_velBC_gimp_damped.ups
```
6. After running for some time, an output file `OneD_pressureBC_from_velBC_gimp_damped.uda.000` will be 
   produced.
7. Run the shell script `extractDataOneD.sh` to extract the relevant data from these MPM runs
   after uncommenting the appropriate lines.  
7.1 If you want to generate a pressure BC from the reactions produced by an impact simulation, copy
    the `RigidReactionForce.dat` file from `OneD_pressureBC_from_velBC_gimp_damped.uda.000` and
    then run the R script `writeReactionLoadCurve.R` to generate a reaction xml file that you can include
    in the `PhysicalBC` section of the input file.

Comparing the 1D and MPM solutions
----------------------------------------
1. You will be able to generate a plot of the exact vs. simulated solution by running the
   R script `plotPulseSim_vs_Exact.R`.   This works only for the step function and needs outputs
   from the velocity BC run and the pressure BC from velocity run.
2. To animate the exact vs. simulated solution try the R script `animatePulseSim_vs_Exact.R`.  This
   is also working only partially.
3. More general animations of simulations can be created by modifying and running the script
   `animatePulseSim.R`.  Animated GIFs of the results are produced by this script.




