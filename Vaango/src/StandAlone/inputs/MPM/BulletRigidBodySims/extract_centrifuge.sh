export WORKING_DIR=.
export TIMEXTRACT_EXE=/home/user/ParSim/Vaango/opt/StandAlone/timeextract
export SELECTPART_EXE=/home/user/ParSim/Vaango/opt/StandAlone/selectpart
export EXTRACTPOS_EXE=/home/user/ParSim/Vaango/opt/StandAlone/extractPos
export EXTRACTPSCA_EXE=/home/user/ParSim/Vaango/opt/StandAlone/extractPscalar
export EXTRACTPVEC_EXE=/home/user/ParSim/Vaango/opt/StandAlone/extractPvec
export EXTRACTPMAT_EXE=/home/user/ParSim/Vaango/opt/StandAlone/extractPmat
export PARTEXTRACT_EXE=/home/user/ParSim/Vaango/opt/StandAlone/partextract
export EXTRACTPOSVEL_EXE=/home/user/ParSim/Vaango/opt/StandAlone/extractPosVelMasVol

# Material 0 = Al. Plate
# Material 1 = Soil
# Material 2 = Hull
# Location of surface : z = 0.352
# Location of hole center: 0, 0.0, 0.301
# Location of top of Al. plate : z = 0.162
# Location of flat-hull plate: z = 0.37
# Area of flat part of v-hull plate: z = 0.0023132 m^2

##./extract_centrifuge_acc.py --x --point 0 0.02 0.162 --mat 0 --uda Centrifuge_Hull_BoulderClay_20g_13ww_midPBC.uda.000

##./extract_centrifuge_acc.py --x --point 0 0.02 0.301 --mat 1 --uda Centrifuge_Hull_BoulderClay_20g_13ww_midPBC.uda.000

#./extract_centrifuge_press.py --x --point 0 0.02 0.301 --dx 0.006 --mat 1 --uda Centrifuge_Hull_BoulderClay_20g_13ww_midPBC.uda.000

#./extract_centrifuge_sigz.py --x --point 0 0.02 0.301 --dx 0.006 --mat 1 --uda Centrifuge_Hull_BoulderClay_20g_13ww_midPBC.uda.000

##./extract_centrifuge_strainrate.py --x --point 0 0.02 0.301 --dx 0.006 --delT 2.39e-7 --mat 1 --uda Centrifuge_Hull_BoulderClay_20g_13ww_midPBC.uda.000

#./extract_centrifuge_bulge.py --mat 1 --line -0.575 0.0 0.352 0.575 0.0 0.352 --tlo 0 --thi 142 --tstep 10 --uda Centrifuge_Hull_BoulderClay_20g_13ww_midPBC.uda.000

##./extract_centrifuge_fext.py --plane -0.5 -0.5 0.37 0.5 -0.5 0.37 -0.5 0.5 0.37 --surf_area 0.007854 --mat 2 --tlo 0 --thi 142 --tstep 10 --uda Centrifuge_Hull_BoulderClay_20g_13ww_midPBC.uda.000

#./extract_centrifuge_mom.py --box -0.5 -0.5 0.36 0.5 0.5 0.47 --mat 2 --tlo 0 --thi 142 --tstep 10 --uda Centrifuge_Hull_BoulderClay_20g_13ww_midPBC.uda.000

#./extract_centrifuge_pos_data.py --box -0.05 -0.05 0.2 0.2  0.2 0.45 --mat 1 --timestep 140 --uda Centrifuge_Hull_BoulderClay_20g_13ww_midPBC.uda.000

#./extract_centrifuge_pos_data.py --box -0.05 -0.05 0.3 0.5  0.5 0.6 --mat 2 --timestep 140 --uda Centrifuge_Hull_BoulderClay_20g_13ww_midPBC.uda.000

./extract_centrifuge_pos_data.py --box 0.25 -0.05 0.2 0.5 0.2 0.45 --mat 1 --timestep 140 --uda Centrifuge_Hull_BoulderClay_20g_13ww_midPBC.uda.000

#./extract_centrifuge_pos_data.py --box -0.05 0.22 0.2 0.5 0.5 0.45 --mat 1 --timestep 140 --uda Centrifuge_Hull_BoulderClay_20g_13ww_midPBC.uda.000

