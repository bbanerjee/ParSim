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

    if (model.attrib['type'] == 'tabular_plasticity_cap'):

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
          if (yield_model.attrib['type'] == 'tabular_cap'):
            filename = yield_model.find('filename').text
            material_dict['yield_filename'] = filename
            R = yield_model.find('cap_ellipticity_ratio').text
            material_dict['R'] = float(R)

        for cap_model in child.iter('internal_variable_model'):
          if (cap_model.attrib['type'] == 'tabular_cap'):
            filename = cap_model.find('filename').text
            material_dict['cap_filename'] = filename

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
# Compute points on the cap
#-----------------------------------------------------------------------------
def bernstein_poly(i, n, t):
  return comb(n, i) * ( t**(n-i) ) * (1 - t)**i

def bezier_curve(xPoints, yPoints, nTimes=1000):
  nPoints = len(xPoints)
  t = np.linspace(0.0, 1.0, nTimes)
  polynomial_array = np.array([ bernstein_poly(i, nPoints-1, t) for i in range(0, nPoints)   ])
  xvals = np.dot(xPoints, polynomial_array)
  yvals = np.dot(yPoints, polynomial_array)
  return xvals, yvals

def computeConvexHull(pbars, sqrtJ2s):

  points = np.ndarray(shape = (len(pbars),2))
  points[:,0] = pbars
  points[:,1] = sqrtJ2s
  if (len(points) < 3):
    hull_points = points
  else:
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
  
def computeEllipseHeight(pbars, sqrtJ2s, p_cap):
  startp = 0;
  endp = 1;
  for i in range(0, len(pbars)):
    if (p_cap > pbars[i]):
      startp = i
      endp = i+1

  if (startp == len(pbars)-1):
    startp = startp - 1
    endp = endp - 1

  t = (p_cap - pbars[startp])/(pbars[endp] - pbars[startp])
  sqrtJ2_cap = (1 - t)*sqrtJ2s[startp] + t*sqrtJ2s[endp]
  #print("t = ", t, "sqrtJ2_cap = ", sqrtJ2_cap)

  return sqrtJ2_cap


def computeCapYieldSurfacePoints(pbars, sqrtJ2s, R, capX):

  # Compute convex hull of tabular data
  hull = computeConvexHull(pbars, sqrtJ2s)
  pbars_hull = list(hull[:,0])
  sqrtJ2s_hull = list(hull[:,1])

  #print("pbar_hull = ", pbars_hull)
  #print("sqrtJ2_hull = ", sqrtJ2s_hull)

  # Find closest point projections of input data to convex hull
  closest = computeClosestPoints(pbars, sqrtJ2s, pbars_hull, sqrtJ2s_hull)
  pbars = list(closest[:,0])
  sqrtJ2s = list(closest[:,1])

  #print("pbar_closest = ", pbars)
  #print("sqrtJ2_closest = ", sqrtJ2s)

  #print("capX = ", capX)

  # Increase density of points
  num_pts = len(pbars)
  
  if (num_pts < 6):
    pbars_new = []
    sqrtJ2s_new = []
    for ii in range(0, num_pts-1):
      x_start = pbars[ii]
      x_end = pbars[ii+1]
      y_start = sqrtJ2s[ii]
      y_end = sqrtJ2s[ii+1]
      dx = (x_end - x_start)/100
      x_vec = np.linspace(x_start, x_end - dx, 100)
      for x in x_vec:
        t = (x - x_start)/(x_end - x_start)
        y = (1 - t)*y_start + t*y_end
        pbars_new.append(x)
        sqrtJ2s_new.append(y)
    pbars_new.append(pbars[num_pts-1])
    sqrtJ2s_new.append(sqrtJ2s[num_pts-1])
    pbars = pbars_new
    sqrtJ2s = sqrtJ2s_new

  #print("pbar_dense = ", pbars)
  #print("sqrtJ2_dense = ", sqrtJ2s)

  # Convert capX to capXbar
  capXbar = -capX

  # Compute kappa
  pbar_min = min(pbars)
  pbar_max = capXbar/3
  kappa = (1 - R)*pbar_min + R*pbar_max

  # Set up the Drucker-Prager cone
  startp = 0
  for i in range(0, len(pbars)):
    if (kappa > pbars[i]):
      startp = i
  if startp == 0:
    startp = 1

  pbar_cone = pbars[0:startp];
  sqrtJ2_cone = sqrtJ2s[0:startp];

  #print("kappa = ", kappa)
  #print("pbar_cone = ", pbar_cone)
  #print("sqrtJ2_cone = ", sqrtJ2_cone)

  # Set up theta
  theta = np.linspace(0, np.pi/2, 10)
  theta = theta[::-1]

  # Set up p_cap
  a = pbar_max - kappa
  pbar_cap = kappa + a*np.cos(theta)

  # Compute sqrtJ2_cap
  sqrtJ2_cap = []
  for i in range(0, len(theta)):
    b = computeEllipseHeight(pbars, sqrtJ2s, pbar_cap[i])
    sqrtJ2_cap.append(b*np.sin(theta[i]))

  #print("pbar_cap = ", pbar_cap)
  #print("sqrtJ2_cap = ", sqrtJ2_cap)

  pbar_yield = pbar_cone + list(pbar_cap)
  sqrtJ2_yield = sqrtJ2_cone + list(sqrtJ2_cap)

  #print("pbar_yield = ", pbar_yield)
  #print("sqrtJ2_yield = ", sqrtJ2_yield)

  # Fit quadratic Bezier curves to the yield polyline
  #pbar_yield, sqrtJ2_yield = bezier_curve(pbar_yield, sqrtJ2_yield, 500)

  return pbar_yield, sqrtJ2_yield
  
#-----------------------------------------------------------------------------
# Plot the yield surface (compression positive)
#-----------------------------------------------------------------------------
def plotPQYieldSurfaceSim(plt, material_dict, yield_table, 
                          capX_list, ev_e_list, ev_p_list, time_list, color_list,
                          pmin, pmax, qmax,
                          compression = 'negative', ):

  # Extract the data from the yield table
  pbars   = yield_table['Pressure']
  sqrtJ2s = yield_table['SqrtJ2']

  for jj, capX in enumerate(capX_list):

    # Compute the yield surface with cap
    R = material_dict["R"]
    pbar_yield, sqrtJ2_yield = computeCapYieldSurfacePoints(pbars, sqrtJ2s, R, capX)

    # Get the plastic strain and set up label
    ev_p_str = str(ev_p_list[jj])
    label_str = 'Plastic vol. strain = ' + ev_p_str

    # Plot
    print('Compression = ', compression)
    if (compression == 'negative'):
      ps = list(map(lambda pbar: -pbar, pbar_yield))
      qs = list(map(lambda sqrtJ2 : np.sqrt(3)*sqrtJ2, sqrtJ2_yield))
    else:
      ps = list(map(lambda pbar: pbar, pbar_yield))
      qs = list(map(lambda sqrtJ2 : np.sqrt(3)*sqrtJ2, sqrtJ2_yield))

    #print("yield surface: pmin = ", pmin, "pmax = ", pmax, "p = ", ps)
    #print("yield surface: qmax = ", qmax, "q = ", qs)
   
    line1 = plt.plot(ps,  qs, '-b',linewidth=1, label = label_str)
    line2 = plt.plot(ps, list(map(lambda q : -q,  qs)), '-b',linewidth=1)  
    plt.setp(line1, color=color_list[jj])
    plt.setp(line2, color=color_list[jj])
    plt.legend(loc=2, prop={'size':8}) 
    #plt.axis('equal')

    axes = plt.gca()
    #axes.set_xlim([1.3*pmin, 1.3*pmax])
    #axes.set_ylim([-1.3*qmax, 1.3*qmax])

  return pmin, qmax
   
