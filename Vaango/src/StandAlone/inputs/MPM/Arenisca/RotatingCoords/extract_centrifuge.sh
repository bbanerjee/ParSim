export WORKING_DIR=.

#./extract_centrifuge_data.py --box -0.2 -0.2 0.217 0.4 0.4 0.25 --mat 0 --timestep 14 --uda Centrifuge_AGR_SimPBC_dense_layer_lores.uda.011
#./extract_centrifuge_data.py --box -0.2 -0.2 0.217 0.4 0.4 0.25 --mat 1 --timestep 14 --uda Centrifuge_AGR_SimPBC_dense_layer_lores.uda.011
#./extract_centrifuge_acc.py --x --point 0 0 0.125 --mat 0 --uda Centrifuge_AGR_SimPBC_dense_layer_lores.uda.009
# ./extract_centrifuge_acc.py --y --point 0 0 0.164 --mat 0 --uda Centrifuge_2D_AGR_SimPBC_dense_layer_lores.uda_loforce.000
#./extract_centrifuge_acc.py --y --point 0 0 0.164 --mat 0 --uda Centrifuge_2D_AGR_SimPBC_dense_layer_lores.uda.001
#./extract_centrifuge_press.py --y --point 0 0 0.164 --dx 0.003 --mat 0 --uda Centrifuge_2D_AGR_SimPBC_dense_layer_lores.uda.001
#./extract_centrifuge_press.py --y --point 0 0 0.164 --dx 0.003 --mat 0 --uda Centrifuge_2D_AGR_SimPBC_dense_layer_lores.uda_loforce.000
#./extract_centrifuge_press.py --x --point 0 0 0.125 --dx 0.003 --mat 0 --uda Centrifuge_AGR_SimPBC_dense_layer_lores.uda.009
#./extract_centrifuge_strainrate.py --y --point 0 0 0.164 --dx 0.003 --mat 0 --uda Centrifuge_2D_AGR_SimPBC_dense_layer_lores.uda.001 --delT 0.001
#./extract_centrifuge_strainrate.py --y --point 0 0 0.164 --dx 0.003 --mat 0 --uda Centrifuge_2D_AGR_SimPBC_dense_layer_lores.uda_loforce.000 --delT 0.001
#./extract_centrifuge_strainrate.py --x --point 0 0 0.125 --dx 0.003 --mat 0 --uda Centrifuge_AGR_SimPBC_dense_layer_lores.uda.009 --delT 0.001

./extract_centrifuge_data.py --box -0.2 -0.2 0.217 0.4 0.4 0.25 --mat 1 --timestep 14 --uda Centrifuge_AGR_SimPBC_dense_layer_very_lores_drained.uda.000
./extract_centrifuge_acc.py --x --point 0 0 0.125 --mat 0 --uda Centrifuge_AGR_SimPBC_dense_layer_very_lores_drained.uda.000
./extract_centrifuge_press.py --x --point 0 0 0.125 --dx 0.006 --mat 0 --uda Centrifuge_AGR_SimPBC_dense_layer_very_lores_drained.uda.000
./extract_centrifuge_strainrate.py --x --point 0 0 0.125 --dx 0.006 --mat 0 --delT 0.001 --uda Centrifuge_AGR_SimPBC_dense_layer_very_lores_drained.uda.000
extract_centrifuge_bulge.py --mat 0 --box -0.575 -0.0 0.20 1.15 0.01 0.01 --tlo 0 --thi 30 --tstep 3 --uda Centrifuge_AGR_SimPBC_dense_layer_very_lores_saturated.uda.000
