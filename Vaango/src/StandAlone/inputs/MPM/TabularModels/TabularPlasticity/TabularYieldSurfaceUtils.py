import os
import xml.etree.ElementTree as ET
import json
import matplotlib.cm as cm
import numpy as np
from scipy.spatial import ConvexHull
from scipy.misc import comb
from shapely.geometry import LineString, Point

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

        for yield_model in child.iter('yield_condition'):
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
    #elem = r'$\mathbf{'+elem+'}$'
    elem = r''+elem+''
    return_string+=elem+'  '
  return return_string[0:-1]
  
#-----------------------------------------------------------------------------
# Compute points on the yield surface
#-----------------------------------------------------------------------------

def computeConvexHull(pbars, sqrtJ2s):

  points = np.ndarray(shape = (len(pbars),2))
  points[:,0] = pbars
  points[:,1] = sqrtJ2s
  hull = ConvexHull(points)

  hull_points = points[hull.vertices]
  hull_points.sort(axis=0)

  return hull_points

def computeClosestPoints(pbars, sqrtJ2s, pbars_hull, sqrtJ2s_hull):

  hull = list(zip(pbars_hull, sqrtJ2s_hull))
  polyline = LineString(hull)

  closest_polyline = []
  for pbar, sqrtJ2 in zip(pbars, sqrtJ2s):
    pt = Point(pbar, sqrtJ2)
    closest = polyline.interpolate(polyline.project(pt))
    closest_polyline.append((closest.x, closest.y))

  return np.array(closest_polyline)

def computeYieldSurfacePoints(pbars, sqrtJ2s):

  # Compute convex hull of tabular data
  if len(pbars) > 3:
    hull = computeConvexHull(pbars, sqrtJ2s)
    pbars_hull = list(hull[:,0])
    sqrtJ2s_hull = list(hull[:,1])
  else:
    pbars_hull = pbars
    sqrtJ2s_hull = sqrtJ2s

  # Find closest point projections of input data to convex hull
  closest = computeClosestPoints(pbars, sqrtJ2s, pbars_hull, sqrtJ2s_hull)
  pbar_yield = list(closest[:,0])
  sqrtJ2_yield = list(closest[:,1])

  return pbar_yield, sqrtJ2_yield
  
#-----------------------------------------------------------------------------
# Plot the yield surface (compression positive)
#-----------------------------------------------------------------------------
def plotPQYieldSurfaceSim(plt, material_dict, yield_table,
                          ev_e_list, ev_p_list, time_list, pmin, pmax, qmax,
                          compression = 'negative', plt_color = 'b'):

  # Extract the data from the yield table
  pbars = yield_table['Pressure']
  sqrtJ2s   = yield_table['SqrtJ2']

  # Compute the convex yield surface
  pbar_yield, sqrtJ2_yield = computeYieldSurfacePoints(pbars, sqrtJ2s)

  # Plot
  #print('Compression = ', compression)
  if (compression == 'negative'):
    ps = list(map(lambda pbar: -pbar, pbar_yield))
    qs = list(map(lambda sqrtJ2 : np.sqrt(3)*sqrtJ2, sqrtJ2_yield))
    #print("yield surface: pmin = ", pmin, "pmax = ", pmax, "p = ", ps)
    #print("yield surface: qmax = ", qmax, "q = ", qs)
 
    line1 = plt.plot(ps,  qs, '-b',linewidth=1)
    line2 = plt.plot(ps, list(map(lambda q : -q,  qs)), '-b',linewidth=1)  
    plt.setp(line1, color=plt_color)
    plt.setp(line2, color=plt_color)
    plt.legend(loc=2, prop={'size':8}) 
    #plt.axis('equal')

    axes = plt.gca()
    #axes.set_xlim([1.3*pmin, 1.3*pmax])
    #axes.set_ylim([-1.3*qmax, 1.3*qmax])
  else:
    ps = list(map(lambda pbar: pbar, pbar_yield))
    qs = list(map(lambda sqrtJ2 : np.sqrt(3)*sqrtJ2, sqrtJ2_yield))
    #print("yield surface: pmin = ", pmin, "pmax = ", pmax, "p = ", ps)
    #print("yield surface: qmax = ", qmax, "q = ", qs)
    line1 = plt.plot(ps,  qs, '-b',linewidth=1)
    line2 = plt.plot(ps, list(map(lambda q : -q,  qs)), '-b',linewidth=1)  
    plt.setp(line1, color=plt_color)
    plt.setp(line2, color=plt_color)
    plt.legend(loc=2, prop={'size':8}) 
    plt.axis('equal')

    axes = plt.gca()
    #axes.set_xlim([1.3*pmin, 1.3*pmax])
    #axes.set_ylim([-1.3*qmax, 1.3*qmax])
 

  return pmin, qmax
   
