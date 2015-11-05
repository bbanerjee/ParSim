#------------------------------------------------------------------------------
#./lineextract -v g.acceleration -startPt -0.5 0 0 -endPt 0.5 0 0 -timestep 200 -uda OneD_tractionBC_gimp.uda.000 > t0.dat
#./lineextract -v g.acceleration -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_tractionBC_gimp.uda.000 > undamped_midres_tractionBC_gacc.dat
#------------------------------------------------------------------------------
#./lineextract -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -timestep 200 -uda OneD_tractionBC_gimp.uda.000 > t0.dat
#./lineextract -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_tractionBC_gimp.uda.000 > undamped_midres_tractionBC_gvel.dat
#------------------------------------------------------------------------------
#./lineextract -v p.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -timestep 200 -uda OneD_tractionBC_gimp.uda.000 > t0.dat
#./lineextract -v p.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_tractionBC_gimp.uda.000 > undamped_midres_tractionBC_pvel.dat
#------------------------------------------------------------------------------
#./lineextract -v p.stress -startPt -0.5 0 0 -endPt 0.5 0 0 -timestep 5 -uda OneD_tractionBC_gimp.uda.000 > t0.dat
#./lineextract -v p.stress -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_tractionBC_gimp.uda.000 > undamped_midres_tractionBC_pstress.dat
#------------------------------------------------------------------------------
#./lineextract -v p.stress -startPt -0.5 0 0 -endPt 0.5 0 0 -timestep 100 -uda OneD_tractionBC_gimp_damped.uda.000 > t0.dat
#./lineextract -v p.stress -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_tractionBC_gimp_damped.uda.000 > damped_midres_tractionBC_pstress.dat
#------------------------------------------------------------------------------
#./lineextract -v p.stress -startPt -0.6 0 0 -endPt 0.6 0 0 -timestep 100 -uda OneD_velocityBC_gimp.uda.000 > t0.dat
#./lineextract -v p.stress -startPt -0.6 0 0 -endPt 0.6 0 0 -uda OneD_velocityBC_gimp.uda.000 > undamped_midres_velBC_pstress.dat
#------------------------------------------------------------------------------
#./lineextract -v p.stress -startPt -0.6 0 0 -endPt 0.6 0 0 -timestep 100 -uda OneD_velocityBC_gimp_damped.uda.000 > t0.dat
#./lineextract -v p.stress -startPt -0.6 0 0 -endPt 0.6 0 0 -uda OneD_velocityBC_gimp_damped.uda.000 > damped_midres_velBC_pstress.dat
#------------------------------------------------------------------------------
#./lineextract -v g.velocity -startPt -0.7 0 0 -endPt 0.7 0 0 -timestep 100 -uda OneD_velocityBC_gimp_damped.uda.000 > t0.dat
#./lineextract -v g.velocity -startPt -0.7 0 0 -endPt 0.7 0 0 -uda OneD_velocityBC_gimp_damped.uda.000 > damped_midres_velBC_gvel.dat
#------------------------------------------------------------------------------
#./lineextract -v p.velocity -startPt -0.6 0 0 -endPt 0.6 0 0 -timestep 100 -uda OneD_velocityBC_gimp_damped.uda.000 > t0.dat
#./lineextract -v p.velocity -startPt -0.6 0 0 -endPt 0.6 0 0 -uda OneD_velocityBC_gimp_damped.uda.000 > damped_midres_velBC_pvel.dat
#------------------------------------------------------------------------------
#./lineextract -v g.velocity -startPt -0.7 0 0 -endPt 0.7 0 0 -timestep 100 -uda OneD_velocityBC_gimp_momform.uda.000 > t0.dat
#./lineextract -v g.velocity -startPt -0.7 0 0 -endPt 0.7 0 0 -uda OneD_velocityBC_gimp_momform.uda.000 > undamped_momform_midres_velBC_gvel.dat
#------------------------------------------------------------------------------
#./lineextract -m 1 -v p.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -timestep 100 -uda OneD_impact_velocityBC_gimp.uda.000 > t0.dat
#./lineextract -m 1 -v p.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_impact_velocityBC_gimp.uda.000 > undamped_midres_impact_velBC_pvel.dat
#------------------------------------------------------------------------------
#./lineextract -v g.velocity -startPt -0.7 0 0 -endPt 0.7 0 0 -timestep 100 -uda OneD_velocityBC_gimp_damped_uintah.uda.000 > t0.dat
#./lineextract -v g.velocity -startPt -0.7 0 0 -endPt 0.7 0 0 -uda OneD_velocityBC_gimp_damped_uintah.uda.000 > damped_uintah_midres_velBC_gvel.dat
#------------------------------------------------------------------------------
#./lineextract -m 1 -v p.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -timestep 100 -uda OneD_impact_velocityBC_gimp_damped.uda.000 > t0.dat
#./lineextract -m 1 -v p.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_impact_velocityBC_gimp_damped.uda.000 > damped_midres_impact_velBC_pvel.dat
#------------------------------------------------------------------------------
./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -timestep 100 -uda OneD_impact_velocityBC_gimp_damped.uda.000 > OneD_impact_velocityBC_gimp_damped.dat
./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_impact_velocityBC_gimp_damped.uda.000 > damped_midres_impact_velBC_gvel.dat
./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -timestep 100 -uda OneD_pressureBC_from_velBC_gimp_damped.uda.000 > OneD_pressureBC_from_velBC_gimp_damped.uda.dat
#------------------------------------------------------------------------------
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -timestep 100 -uda OneD_impact_velocityBC_hat_gimp_damped.uda.000 > OneD_impact_velocityBC_hat_gimp_damped.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -timestep 100 -uda OneD_pressureBC_from_velBC_hat_gimp_damped.uda.000 > OneD_pressureBC_from_velBC_hat_gimp_damped.dat
#------------------------------------------------------------------------------
# Extract all timesteps
#------------------------------------------------------------------------------
#./lineextract -m 1 -v p.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_impact_velocityBC_gimp_damped.uda.000 > damped_midres_SquareVelBC_pvel.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_impact_velocityBC_gimp_damped.uda.000 > damped_midres_SquareVelBC.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_from_velBC_gimp_damped.uda.000 > damped_midres_SquarePressFromVelBC.dat
#------------------------------------------------------------------------------
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_impact_velocityBC_hat_gimp_damped.uda.001 > damped_midres_HatVelBC.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_from_velBC_hat_gimp_damped.uda.000 > damped_midres_PressFromVelBC.dat
#------------------------------------------------------------------------------
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_from_velBC_hat_gimp_undamped_arenisca_elastic_lores.uda.000 > undamped_arenisca_lores_PressFromVelBC.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_from_velBC_hat_gimp_damped_arenisca_elastic_lores.uda.000 > damped_arenisca_lores_PressFromVelBC.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_from_velBC_hat_gimp_undamped_arenisca_elastic.uda.000 > undamped_arenisca_midres_PressFromVelBC.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_from_velBC_hat_gimp_damped_arenisca_elastic.uda.000 > damped_arenisca_midres_PressFromVelBC.dat
#------------------------------------------------------------------------------
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_explosion_gimp_undamped.uda.000 > undamped_midres_explosionBC.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_explosion_gimp_damped.uda.000 > damped_midres_explosionBC.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_explosion_gimp_undamped_lores.uda.000 > undamped_lores_explosionBC.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_explosion_gimp_damped_lores.uda.000 > damped_lores_explosionBC.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_explosion_gimp_damped_hires.uda.000 > damped_hires_explosionBC.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_explosion_gimp_damped_lores_cfl.uda.000 > damped_lores_cfl_explosionBC.dat
#./lineextract -m 1 -v g.velocity -startPt -0.5 0 0 -endPt 0.5 0 0 -uda OneD_pressureBC_explosion_cpdi_damped_lores.uda.000 > damped_lores_cpdi_explosionBC.dat
#------------------------------------------------------------------------------
