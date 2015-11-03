#
# Shell script to link input files and run Aldridge solution code
#
rm fort.1 fort.2 fort.3
ln -s ../input_param_file_hires.in fort.1
ln -s ../CentrifugeLoadCurveLow.xml.resampled fort.2
./explod_for.exe 
cp fort.3 output_data_hires.out
#
rm fort.1 fort.2 fort.3
ln -s ../input_param_file_lores.in fort.1
ln -s ../CentrifugeLoadCurveLow.xml.resampled fort.2
./explod_for.exe 
cp fort.3 output_data_lores.out
