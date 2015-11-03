export WORKING_DIR=.
export TIMEXTRACT_EXE=/home/banerjee/ParSim/Vaango/opt/StandAlone/timeextract
export SELECTPART_EXE=/home/banerjee/ParSim/Vaango/opt/StandAlone/selectpart
export EXTRACTPOS_EXE=/home/banerjee/ParSim/Vaango/opt/StandAlone/extractPos
export PARTEXTRACT_EXE=/home/banerjee/ParSim/Vaango/opt/StandAlone/partextract
export STRESSEXTRACT_EXE=/home/banerjee/ParSim/Vaango/opt/StandAlone/extractS

#./extract_fullspace_acc.py --point 0 0 0.60 --theta 0 --phi 0 --mat 0 --uda Spherical_source_verification_gimp.uda.000

#./extract_fullspace_press.py --x --point 0 0 0.60 --dx 0.006 --mat 0 --uda Spherical_source_verification_gimp.uda.000

./extract_fullspace_partdata.py --point 0 0 0.60 --theta 0 --phi 0 --dx 0.006 --mat 0 --uda Spherical_source_verification_gimp.uda.000

#./extract_fullspace_acc.py --point 0 0 0.60 --theta 0 --phi 0 --mat 0 --uda Spherical_source_verification_gimp_quarter.uda.000

#./extract_fullspace_press.py --x --point 0 0 0.60 --dx 0.006 --mat 0 --uda Spherical_source_verification_gimp_quarter.uda.000

#./extract_fullspace_partdata.py --point 0 0 0.60 --theta 0 --phi 0 --dx 0.006 --mat 0 --uda Spherical_source_verification_gimp_quarter.uda.000

#./extract_fullspace_acc.py --point 0 0 0.60 --theta 0 --phi 0 --mat 0 --uda Spherical_source_verification_gimp_quarter_hires.uda.000

#./extract_fullspace_partdata.py --point 0 0 0.60 --theta 0 --phi 0 --dx 0.006 --mat 0 --uda Spherical_source_verification_gimp_quarter_hires.uda.000

#./extract_fullspace_acc.py --point 0 0 0.60 --theta 0 --phi 0 --mat 0 --uda Spherical_source_verification_gimp_quarter_hires_cpdi.uda.000

#./extract_fullspace_partdata.py --point 0 0 0.60 --theta 0 --phi 0 --dx 0.006 --mat 0 --uda Spherical_source_verification_gimp_quarter_hires_cpdi.uda.000

#./extract_fullspace_acc.py --point 0 0 0.60 --theta 0 --phi 0 --mat 0 --uda Spherical_source_verification_gimp_quarter_hires_cpdi_nocbdi.uda.000

#./extract_fullspace_partdata.py --point 0 0 0.60 --theta 0 --phi 0 --dx 0.006 --mat 0 --uda Spherical_source_verification_gimp_quarter_hires_cpdi_nocbdi.uda.000

#./extract_fullspace_acc.py --point 0 0 0.60 --theta 0 --phi 0 --mat 0 --uda Spherical_source_verification_gimp_quarter_hires_cpdi_nocbdi_nodamp.uda.000

#./extract_fullspace_partdata.py --point 0 0 0.60 --theta 0 --phi 0 --dx 0.006 --mat 0 --uda Spherical_source_verification_gimp_quarter_hires_cpdi_nocbdi_nodamp.uda.000
