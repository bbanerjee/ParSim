from TabularTest_01_UniaxialStrainRotateJ2Lin import *

def post_proc(test, uda_path, save_path, POST_PROCESS_LIST):
  test_name = os.path.split(test)[1]
  if test_name in POST_PROCESS_LIST:
    if test_name == 'TabularTest_01_UniaxialStrainRotateJ2Lin.ups':
      uniaxialStrainRotateJ2Lin(uda_path,save_path)
  else:
    print('\nERROR: test: ',test,'\n\tNot on post processing list.\n')

