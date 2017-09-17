import os
import xml.etree.ElementTree as ET
import json
import matplotlib.cm as cm
import numpy as np

#-----------------------------------------------------------------------------
# Read the input xml file and save relevant data to a dictionary
#-----------------------------------------------------------------------------
def get_yield_surface_data(uda_path, PRINTOUT=False):

  # Parse xml
  try:
    ups_file = os.path.abspath(uda_path)+'/input.xml.orig'
    tree = ET.parse(ups_file)
  except:
    ups_file = os.path.abspath(uda_path)+'/input.xml'
    tree = ET.parse(ups_file)

  # Read the plasticity parameters
  material_dict = {}
  root = tree.getroot()
  for model in root.iter('constitutive_model'):

    if (model.attrib['type'] == 'tabular_plasticity'):

      for child in model:

        for elastic_model in child.iter('elastic_moduli_model'):
          if (elastic_model.attrib['type'] == 'tabular'):
            filename = elastic_model.find('filename').text
            G0 = elastic_model.find('G0').text
            nu = elastic_model.find('nu').text
            material_dict['elastic_filename'] = filename
            material_dict['G0'] = float(G0)
            material_dict['nu'] = float(nu)

        for yield_model in child.iter('plastic_yield_condition'):
          if (yield_model.attrib['type'] == 'tabular'):
            filename = yield_model.find('filename').text
            material_dict['yield_filename'] = filename

        scaling = model.find('yield_surface_radius_scaling_factor').text 
        subcycling = model.find('subcycling_characteristic_number').text
        material_dict['yield scale fac'] = float(scaling)
        material_dict['subcycling char num']  = int(subcycling)

  if PRINTOUT:
    print('--Material Specification--')
    for key in material_dict:
      print(key,':',material_dict[key])

  material_dict = convertToBoldFont(material_dict, PRINTOUT)

  return material_dict
      
#-----------------------------------------------------------------------------
# Convert material dictionary data to bold font and add string to material_dict
#-----------------------------------------------------------------------------
def convertToBoldFont(material_dict, PRINTOUT=False):

  tmp_string = r'$\mathbf{\underline{Material\phantom{1}Properties:}}$'+'\n'
  key_list = list(material_dict.keys())
  key_list.sort()
  for key in key_list:
    if '_' in key:
      tmp = key.split('_')
      tmp = str_to_mathbf(tmp[0]+'_'+'{'+tmp[1]+'}')
      if (isinstance(material_dict[key], str)):
        tmp_string += tmp + str_to_mathbf(' = ')+str_to_mathbf(material_dict[key])+'\n'
      else:
        tmp_string += tmp+str_to_mathbf(' = ')+str_to_mathbf(format(material_dict[key],'1.3e'))+'\n'
    else:
      tmp = key
      if key == 'subcycling char num':
        tmp_string += str_to_mathbf(tmp+' = '+format(material_dict[key],'4.1f'))+'\n'
      else:
        if (isinstance(material_dict[key], str)):
          tmp_string += str_to_mathbf(tmp+' = '+ material_dict[key])+'\n'
        else:
          tmp_string += str_to_mathbf(tmp+' = '+format(material_dict[key],'1.3e'))+'\n'

  material_dict['material string'] = tmp_string[0:-1]
  if PRINTOUT:
    print(tmp_string)

  return material_dict

#-----------------------------------------------------------------------------
# Make the text bold
# Only works with single spaces no leading space
#-----------------------------------------------------------------------------
def str_to_mathbf(string):
  string = string.split()
  return_string = ''
  for elem in string:
    elem = r'$\mathbf{'+elem+'}$'
    return_string+=elem+'  '
  return return_string[0:-1]
  
#-----------------------------------------------------------------------------
# Plot the yield surface
#-----------------------------------------------------------------------------
def plotPQYieldSurfaceSim(plt, material_dict, yield_table, ev_e_list, ev_p_list, time_list):

  # Extract the data from the yield table
  pressures = yield_table['Pressure']
  sqrtJ2s   = yield_table['SqrtJ2']

  # Convert into p and q (Pa)
  ps = list(map(lambda pbar: -pbar, pressures))
  qs = list(map(lambda sqrtJ2 : np.sqrt(3)*sqrtJ2, sqrtJ2s))
 
  # Choose the Paired colormap
  plt_color = cm.Paired(float(1)/len(ev_p_list))

  # Plot
  line1 = plt.plot(ps,  qs, '-b',linewidth=1)
  line2 = plt.plot(ps, list(map(lambda q : -q,  qs)), '-b',linewidth=1)  
  plt.setp(line1, color=plt_color)
  plt.setp(line2, color=plt_color)
  plt.legend(loc=2, prop={'size':8}) 

  # Find min and max
  pmin = min(ps)
  pmax = max(ps)
  qmax = max(qs)

  axes = plt.gca()
  axes.set_xlim([1.0e-6*pmin, 1.0e-6*pmax])
  #axes.set_xlim([1.3*pmin, 1.3*pmax])
  axes.set_ylim([-1.3*qmax, 1.3*qmax])
  return pmin, qmax
   
