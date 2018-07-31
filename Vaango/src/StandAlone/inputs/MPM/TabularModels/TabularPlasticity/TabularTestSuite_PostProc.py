from TabularTest_01_UniaxialStrainRotateJ2Lin import *
from TabularTest_01a_UniaxialStressJ2Lin import *
from TabularTest_01b_UniaxialStrainLoadUnloadJ2Lin import *
from TabularTest_02_UniaxialStrainLoadUnloadJ2NonLin import *
from TabularTest_03_UniaxialStrainRotateDPLin import *
from TabularTest_03a_UniaxialStrainLoadUnloadDPLin import *
from TabularTest_04_UniaxialStrainLoadUnloadDPNonLin import *
from TabularTest_05_UniaxialStrainLoadUnloadNonLinDPNonLin import *
from TabularTest_06_HydrostaticCompression import *
from TabularTest_06_HydrostaticCompressionNN import *
from TabularTest_07_HydrostaticLoadUnload import *
from TabularTest_07_HydrostaticLoadUnloadNN import *
from TabularTest_08_UniaxialStrainCompresson import *
from TabularTest_08_UniaxialStrainCompressonNN import *
from TabularTest_09_UniaxialStrainTension import *
from TabularTest_09_UniaxialStrainTensionNN import *
from TabularTest_10_UniaxialStrainRotate import *
from TabularTest_10_UniaxialStrainRotateNN import *
from TabularTest_11_TriaxialStrainTension import *
from TabularTest_11_TriaxialStrainTensionNN import *
from TabularTest_12_UniaxialStrainLoadUnload import *
from TabularTest_12_UniaxialStrainLoadUnloadNN import *
from TabularTest_13_MultiaxialStrainLoadUnload import *
from TabularTest_13_MultiaxialStrainLoadUnloadNN import *

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
    elif test_name == 'TabularTest_04_UniaxialStrainLoadUnloadDPNonLin.ups':
      uniaxialStrainLoadUnloadDPNonLin(uda_path,save_path)
    elif test_name == 'TabularTest_05_UniaxialStrainLoadUnloadNonLinDPNonLin.ups':
      uniaxialStrainLoadUnloadNonLinDPNonLin(uda_path,save_path)
    elif test_name == 'TabularTest_06_HydrostaticCompression.ups':
      hydrostaticCompression(uda_path, save_path)
    elif test_name == 'TabularTest_06_HydrostaticCompressionNN.ups':
      hydrostaticCompressionNN(uda_path, save_path)
    elif test_name == 'TabularTest_07_HydrostaticLoadUnload.ups':
      hydrostaticLoadUnload(uda_path, save_path)
    elif test_name == 'TabularTest_07_HydrostaticLoadUnloadNN.ups':
      hydrostaticLoadUnloadNN(uda_path, save_path)
    elif test_name == 'TabularTest_08_UniaxialStrainCompresson.ups':
      uniaxialStrainCompression(uda_path, save_path)
    elif test_name == 'TabularTest_08_UniaxialStrainCompressonNN.ups':
      uniaxialStrainCompressionNN(uda_path, save_path)
    elif test_name == 'TabularTest_09_UniaxialStrainTension.ups':
      uniaxialStrainTension(uda_path, save_path)
    elif test_name == 'TabularTest_09_UniaxialStrainTensionNN.ups':
      uniaxialStrainTensionNN(uda_path, save_path)
    elif test_name == 'TabularTest_10_UniaxialStrainRotate.ups':
      uniaxialStrainTensionNN(uda_path, save_path)
    elif test_name == 'TabularTest_10_UniaxialStrainRotateNN.ups':
      uniaxialStrainRotateNN(uda_path, save_path)
    elif test_name == 'TabularTest_11_TriaxialStrainTension.ups':
      triaxialStrainTension(uda_path, save_path)
    elif test_name == 'TabularTest_11_TriaxialStrainTensionNN.ups':
      triaxialStrainTensionNN(uda_path, save_path)
    elif test_name == 'TabularTest_12_UniaxialStrainLoadUnload.ups':
      uniaxialStrainLoadUnload(uda_path, save_path)
    elif test_name == 'TabularTest_12_UniaxialStrainLoadUnloadNN.ups':
      uniaxialStrainLoadUnloadNN(uda_path, save_path)
    elif test_name == 'TabularTest_13_MultiaxialStrainLoadUnload.ups':
      multiaxialStrainLoadUnload(uda_path, save_path)
    elif test_name == 'TabularTest_13_MultiaxialStrainLoadUnloadNN.ups':
      multiaxialStrainLoadUnloadNN(uda_path, save_path)
    else:
      print('\nERROR: test: ',test,'\n\tNot on post processing list.\n')
  else:
    print('\nERROR: test: ',test,'\n\tNot on post processing list.\n')

