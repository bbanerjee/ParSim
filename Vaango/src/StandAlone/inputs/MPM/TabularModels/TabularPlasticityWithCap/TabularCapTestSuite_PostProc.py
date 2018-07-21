from TabularCapTest_01_HydrostaticCompression import *
from TabularCapTest_01_HydrostaticCompressionNN import *
from TabularCapTest_02_HydrostaticLoadUnload import *
from TabularCapTest_02_HydrostaticLoadUnloadNN import *
from TabularCapTest_03_UniaxialStrainCompresson import *
from TabularCapTest_03_UniaxialStrainCompressonNN import *
from TabularCapTest_04_UniaxialStrainTension import *
from TabularCapTest_04_UniaxialStrainTensionNN import *
from TabularCapTest_05_UniaxialStrainRotate import *
from TabularCapTest_05_UniaxialStrainRotateNN import *
from TabularCapTest_06_TriaxialStrainTension import *
from TabularCapTest_06_TriaxialStrainTensionNN import *
from TabularCapTest_07_UniaxialStrainLoadUnload import *
from TabularCapTest_07_UniaxialStrainLoadUnloadNN import *
from TabularCapTest_08_MultiaxialStrainLoadUnload import *
from TabularCapTest_08_MultiaxialStrainLoadUnloadNN import *

def post_proc(test, uda_path, save_path, POST_PROCESS_LIST):
  test_name = os.path.split(test)[1]
  if test_name in POST_PROCESS_LIST:
    if test_name == 'TabularCapTest_01_HydrostaticCompression.ups':
      hydrostaticCompression(uda_path, save_path)
    elif test_name == 'TabularCapTest_01_HydrostaticCompressionNN.ups':
      hydrostaticCompressionNN(uda_path, save_path)
    elif test_name == 'TabularCapTest_02_HydrostaticLoadUnload.ups':
      hydrostaticLoadUnload(uda_path, save_path)
    elif test_name == 'TabularCapTest_02_HydrostaticLoadUnloadNN.ups':
      hydrostaticLoadUnloadNN(uda_path, save_path)
    elif test_name == 'TabularCapTest_03_UniaxialStrainCompresson.ups':
      uniaxialStrainCompression(uda_path, save_path)
    elif test_name == 'TabularCapTest_03_UniaxialStrainCompressonNN.ups':
      uniaxialStrainCompressionNN(uda_path, save_path)
    elif test_name == 'TabularCapTest_04_UniaxialStrainTension.ups':
      uniaxialStrainTension(uda_path, save_path)
    elif test_name == 'TabularCapTest_04_UniaxialStrainTensionNN.ups':
      uniaxialStrainTensionNN(uda_path, save_path)
    elif test_name == 'TabularCapTest_05_UniaxialStrainRotate.ups':
      uniaxialStrainTensionNN(uda_path, save_path)
    elif test_name == 'TabularCapTest_05_UniaxialStrainRotateNN.ups':
      uniaxialStrainRotateNN(uda_path, save_path)
    elif test_name == 'TabularCapTest_06_TriaxialStrainTension.ups':
      triaxialStrainTension(uda_path, save_path)
    elif test_name == 'TabularCapTest_06_TriaxialStrainTensionNN.ups':
      triaxialStrainTensionNN(uda_path, save_path)
    elif test_name == 'TabularCapTest_07_UniaxialStrainLoadUnload.ups':
      uniaxialStrainLoadUnload(uda_path, save_path)
    elif test_name == 'TabularCapTest_07_UniaxialStrainLoadUnloadNN.ups':
      uniaxialStrainLoadUnloadNN(uda_path, save_path)
    elif test_name == 'TabularCapTest_08_MultiaxialStrainLoadUnload.ups':
      multiaxialStrainLoadUnload(uda_path, save_path)
    elif test_name == 'TabularCapTest_08_MultiaxialStrainLoadUnloadNN.ups':
      multiaxialStrainLoadUnloadNN(uda_path, save_path)
    else:
      print('\nERROR: test: ',test,'\n\tNot on post processing list.\n')
  else:
    print('\nERROR: test: ',test,'\n\tNot on post processing list.\n')

