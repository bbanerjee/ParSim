from J2Test_01_UniaxialStrainRotateJ2Lin import *
from J2Test_01a_UniaxialStressJ2Lin import *
from J2Test_01b_UniaxialStrainLoadUnloadJ2Lin import *
from J2Test_02_UniaxialStrainLoadUnloadJ2NonLin import *
from J2Test_03_UniaxialStrainRotateDPLin import *
from J2Test_03a_UniaxialStrainLoadUnloadDPLin import *
from J2Test_04_UniaxialStrainLoadUnloadDPNonLin import *
from J2Test_05_UniaxialStrainLoadUnloadNonLinDPNonLin import *
from J2Test_06_HydrostaticCompression import *
from J2Test_07_HydrostaticLoadUnload import *
from J2Test_08_UniaxialStrainCompresson import *
from J2Test_09_UniaxialStrainTension import *
from J2Test_10_UniaxialStrainRotate import *
from J2Test_11_TriaxialStrainTension import *
from J2Test_12_UniaxialStrainLoadUnload import *
from J2Test_13_MultiaxialStrainLoadUnload import *

def post_proc(test, uda_path, save_path, POST_PROCESS_LIST):
  test_name = os.path.split(test)[1]
  if test_name in POST_PROCESS_LIST:
    if test_name == 'J2Test_01_UniaxialStrainRotateJ2Lin.ups':
      uniaxialStrainRotateJ2Lin(uda_path,save_path)
    elif test_name == 'J2Test_01a_UniaxialStressJ2Lin.ups':
      uniaxialStressJ2Lin(uda_path,save_path)
    elif test_name == 'J2Test_01b_UniaxialStrainLoadUnloadJ2Lin.ups':
      uniaxialStrainLoadUnloadJ2Lin(uda_path,save_path)
    elif test_name == 'J2Test_02_UniaxialStrainLoadUnloadJ2NonLin.ups':
      uniaxialStrainLoadUnloadJ2NonLin(uda_path,save_path)
    elif test_name == 'J2Test_03_UniaxialStrainRotateDPLin.ups':
      uniaxialStrainRotateDPLin(uda_path,save_path)
    elif test_name == 'J2Test_03a_UniaxialStrainLoadUnloadDPLin.ups':
      uniaxialStrainLoadUnloadDPLin(uda_path,save_path)
    elif test_name == 'J2Test_04_UniaxialStrainLoadUnloadDPNonLin.ups':
      uniaxialStrainLoadUnloadDPNonLin(uda_path,save_path)
    elif test_name == 'J2Test_05_UniaxialStrainLoadUnloadNonLinDPNonLin.ups':
      uniaxialStrainLoadUnloadNonLinDPNonLin(uda_path,save_path)
    elif test_name == 'J2Test_06_HydrostaticCompression.ups':
      hydrostaticCompression(uda_path, save_path)
    elif test_name == 'J2Test_07_HydrostaticLoadUnload.ups':
      hydrostaticLoadUnload(uda_path, save_path)
    elif test_name == 'J2Test_08_UniaxialStrainCompresson.ups':
      uniaxialStrainCompression(uda_path, save_path)
    elif test_name == 'J2Test_09_UniaxialStrainTension.ups':
      uniaxialStrainTension(uda_path, save_path)
    elif test_name == 'J2Test_11_TriaxialStrainTension.ups':
      triaxialStrainTension(uda_path, save_path)
    elif test_name == 'J2Test_12_UniaxialStrainLoadUnload.ups':
      uniaxialStrainLoadUnload(uda_path, save_path)
    elif test_name == 'J2Test_13_MultiaxialStrainLoadUnload.ups':
      multiaxialStrainLoadUnload(uda_path, save_path)
    else:
      print('\nERROR: test: ',test,'\n\tNot on post processing list.\n')
  else:
    print('\nERROR: test: ',test,'\n\tNot on post processing list.\n')

