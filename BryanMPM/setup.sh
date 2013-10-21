#bin/bash
#If not using python 2.7, change include directory to correct location
cd ./src/shapes
cython gimp2_c.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 -o gimp2_c.so gimp2_c.c
cython quad2_c.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 -o quad2_c.so quad2_c.c
cython linear2_c.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 -o linear2_c.so linear2_c.c
cd ..
cython materialmodel2d_c.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 -o materialmodel2d_c.so materialmodel2d_c.c
cython mpmutils_c.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 -o mpmutils_c.so mpmutils_c.c
cd ../examples
mkdir test_data
cd ..
