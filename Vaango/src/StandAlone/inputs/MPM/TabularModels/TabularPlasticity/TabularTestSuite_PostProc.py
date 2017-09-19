from TabularTest_01_UniaxialStrainRotateJ2Lin import *
from TabularTest_01a_UniaxialStressJ2Lin import *
from TabularTest_01b_UniaxialStrainLoadUnloadJ2Lin import *
from TabularTest_02_UniaxialStrainLoadUnloadJ2NonLin import *
from TabularTest_03_UniaxialStrainRotateDPLin import *
from TabularTest_03a_UniaxialStrainLoadUnloadDPLin import *

def post_proc(test, uda_path, save_path, POST_PROCESS_LIST):
  test_name = os.path.split(test)[1]
  if test_name in POST_PROCESS_LIST:
    if test_name == 'TabularTest_01_UniaxialStrainRotateJ2Lin.ups':
      uniaxialStrainRotateJ2Lin(uda_path,save_path)
    elif test_name == 'TabularTest_01a_UniaxialStressJ2Lin.ups':
      uniaxialStressJ2Lin(uda_path,save_path)
    elif test_name == 'TabularTest_01b_UniaxialStrainLoadUnloadJ2Lin.ups':
      uniaxialStrainLoadUnloadJ2Lin(uda_path,save_path)
    elif test_name == 'TabularTest_02_UniaxialStrainLoadUnloadJ2NonLin.ups':
      uniaxialStrainLoadUnloadJ2NonLin(uda_path,save_path)
    elif test_name == 'TabularTest_03_UniaxialStrainRotateDPLin.ups':
      uniaxialStrainRotateDPLin(uda_path,save_path)
    elif test_name == 'TabularTest_03a_UniaxialStrainLoadUnloadDPLin.ups':
      uniaxialStrainLoadUnloadDPLin(uda_path,save_path)
    else:
      print('\nERROR: test: ',test,'\n\tNot on post processing list.\n')
  else:
    print('\nERROR: test: ',test,'\n\tNot on post processing list.\n')

