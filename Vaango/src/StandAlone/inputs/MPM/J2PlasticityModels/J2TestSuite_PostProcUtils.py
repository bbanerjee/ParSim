#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import tempfile
import numpy as np
import subprocess as sub_proc

#Plotting stuff below
from matplotlib import rc
from matplotlib import ticker
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from J2YieldSurfaceUtils import *

SHOW_ON_MAKE = False

# Get the partextact executable path from the environment
partextract_exe = os.path.abspath(os.environ['PARTEXTRACT_EXE'])

#Useful constants
sqrtThree = np.sqrt(3.0)
twoThirds = 2.0 / 3.0
threeHalves = 3.0 / 2.0

#Set matplotlib defaults to desired values
#Set the legend to best fit
fontSize = 16

markers = None
plt.rcParams['legend.loc'] = 'best'
#Set font size
plt.rcParams['mathtext.it'] = 'serif:bold'
plt.rcParams['mathtext.rm'] = 'serif:bold'
plt.rcParams['mathtext.sf'] = 'serif:bold'
plt.rcParams['font.size'] = fontSize
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelsize'] = 'medium'
#plt.rcParams['axes.labelweight']='bold'
plt.rcParams['legend.fontsize'] = 'medium'
#Set linewidth
lineWidth = 2
plt.rcParams['lines.linewidth'] = lineWidth
#Set markersize
plt.rcParams['lines.markersize'] = 8
#Set padding for tick labels and size
plt.rcParams['xtick.major.pad'] = 12
plt.rcParams['ytick.major.pad'] = 8
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3

#resolution
plt.rcParams['figure.dpi'] = 120

font = {'family': 'serif', 'weight': 'bold', 'size': fontSize}
rc('font', **font)
rc('text', usetex=True)


def savePNG(name, size='1920x1080'):
  res = float(plt.rcParams['figure.dpi'])
  #Add Check for file already existing as name.png
  if size == '640x480':
    size = [640 / res, 480 / res]
  if size == '1080x768':
    size = [1080 / res, 768 / res]
  if size == '1152x768':
    size = [1152 / res, 768 / res]
  if size == '1280x854':
    size = [1280 / res, 854 / res]
  if size == '1280x960':
    size = [1280 / res, 960 / res]
  if size == '1920x1080':
    size = [1920 / res, 1080 / res]
  #set the figure size for saving
  #plt.gcf().set_size_inches(size[0],size[1])
  #save at speciified resolution
  #plt.savefig(name+'.png', bbox_inches=0, dpi=plt.rcParams['figure.dpi'])

  #size='960x960'
  plt.gcf().set_size_inches(size[0], size[1])
  plt.savefig(name + '.pdf', bbox_inches=0, dpi=plt.rcParams['figure.dpi'])


def savePDF(name, size='1920x1080'):
  res = float(plt.rcParams['figure.dpi'])
  #Add Check for file already existing as name.png
  if size == '640x480':
    size = [640 / res, 480 / res]
  if size == '1080x768':
    size = [1080 / res, 768 / res]
  if size == '1152x768':
    size = [1152 / res, 768 / res]
  if size == '1280x854':
    size = [1280 / res, 854 / res]
  if size == '1280x960':
    size = [1280 / res, 960 / res]
  if size == '1920x1080':
    size = [1920 / res, 1080 / res]
  #set the figure size for saving
  plt.gcf().set_size_inches(size[0], size[1])
  #save at speciified resolution
  plt.savefig(name + '.pdf', bbox_inches=0, dpi=plt.rcParams['figure.dpi'])

#-----------------------------------------------------------------------------
# Make the text bold
# Only works with single spaces no leading space
#-----------------------------------------------------------------------------
def str_to_mathbf(string):
  string = string.split()
  return_string = ''
  for elem in string:
    #elem = r'$\mathbf{'+elem+'}$'
    elem = r'' + elem + ''
    return_string += elem + '  '
  return return_string[0:-1]


def sign(x, y):
  if y >= 0:
    return abs(x)
  else:
    return -abs(x)


def sigma_iso(sigma):
  return (np.trace(sigma) / 3.0) * np.eye(3)


def sigma_dev(sigma):
  return sigma - sigma_iso(sigma)


def sigma_I1(sigma):
  #print(sigma)
  return sigma.trace()


def sigma_J2(sigma):
  return 0.5 * np.dot(sigma_dev(sigma), sigma_dev(sigma)).trace()


def sigma_J3(sigma):
  return (1 / 3.0) * np.dot(np.dot(sigma_dev(sigma), sigma_dev(sigma)),
                            sigma_dev(sigma)).trace()


def sigma_mag(sigma):
  #Returns the magnitude of a second-rank tensor
  #return np.linalg.norm(sigma)
  return np.sqrt(DblDot(sigma, sigma))


def DblDot(x, y):  #Returns the double inner product of two second-rank tensors
  val = 0
  for i in range(0, 3):
    for j in range(0, 3):
      val = val + (x[i][j] * y[i][j])
  return val


def sigma_tau(sigma):
  #return sign(np.sqrt(sigma_J2(sigma)),sigma_J3(sigma))
  return sign(np.sqrt(sigma_J2(sigma)), sigma_J3(sigma))

def tensor_det(tensor):
  return np.linalg.det(tensor)


def get_ps_and_qs(sigmas):
  ps = []
  qs = []
  for sigma in sigmas:
    qs.append(sign(sqrtThree * np.sqrt(sigma_J2(sigma)), sigma_J3(sigma)))
    ps.append(sigma_I1(sigma) / 3.0)
  return ps, qs


#---------------------------------------------------------------------------------
# Read the internal state variables from the UDA file for a set of times
# using linear interpolation
#---------------------------------------------------------------------------------
def getInternalVariables(uda_path, analytical_times, matID=0):

  # Get the internal variables
  eps_p_eq, epsdot_p_eq, backstress, porosity, damage, temp, vol, times = getAllInternalVariables(
      uda_path, matID)

  # Create a clean list containing the data as functions of time
  eps_p_list = []
  epsdot_p_list = []
  backstress_list = []
  porosity_list = []
  damage_list = []
  temp_list = []
  vol_list = []
  time_list = []

  #print(analytical_times)
  #print(times)
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  #print("Analytical times = ", an_times_add)
  time_list.append(times[0])
  eps_p_list.append(eps_p_eq[0])
  epsdot_p_list.append(epsdot_p_eq[0])
  backstress_list.append(backstress[0])
  porosity_list.append(porosity[0])
  damage_list.append(damage[0])
  temp_list.append(temp[0])
  vol_list.append(vol[0])
  for ii in range(1, len(times)):
    tt_0 = times[ii - 1]
    tt_1 = times[ii]
    for jj, ta in enumerate(an_times_add):
      ss = (ta - tt_0) / (tt_1 - tt_0)
      #print("ii = " , ii, " tt_0 = " , tt_0, " tt_1 = ", tt_1, " jj = " , jj, " ta = " , ta, " ss = ", ss )
      if (ss > 0.0 and ss < 1.0):
        #print("ii = " , ii, " tt_0 = " , tt_0, " tt_1 = ", tt_1, " jj = " , jj, " ta = " , ta )
        time_list.append(ta)
        eps_p_list.append((1 - ss) * eps_p_eq[ii - 1] + ss * eps_p_eq[ii])
        epsdot_p_list.append((1 - ss) * epsdot_p_eq[ii - 1] +
                             ss * epsdot_p_eq[ii])
        backstress_list.append((1 - ss) * backstress[ii - 1] +
                               ss * backstress[ii])
        porosity_list.append((1 - ss) * porosity[ii - 1] + ss * porosity[ii])
        damage_list.append((1 - ss) * damage[ii - 1] + ss * damage[ii])
        temp_list.append((1 - ss) * temp[ii - 1] + ss * temp[ii])
        vol_list.append((1 - ss) * vol[ii - 1] + ss * vol[ii])
    #    break
    #else:
    #  continue
  time_list.append(times[len(times) - 1])
  eps_p_list.append(eps_p_eq[len(times) - 1])
  epsdot_p_list.append(epsdot_p_eq[len(times) - 1])
  backstress_list.append(backstress[len(times) - 1])
  porosity_list.append(porosity[len(times) - 1])
  damage_list.append(damage[len(times) - 1])
  temp_list.append(temp[len(times) - 1])
  vol_list.append(vol[len(times) - 1])

  #print("time = ", time_list)
  #print("eps_p_eq = ", eps_p_list)
  #print("epsdot_p_eq = ", epsdot_p_list)
  #print("backstress = ", backstress_list)
  #print("porosity = ", porosity_list)
  #print("damage = ", damage_list)
  #print("temperature = ", temp_list)
  #print("volume = ", vol_list)

  return eps_p_list, epsdot_p_list, backstress_list, porosity_list, damage_list, temp_list, vol_list, time_list


#---------------------------------------------------------------------------------
# Read the internal state variables from the UDA file for a set of times
# using closest point interpolation
#---------------------------------------------------------------------------------
def find_nearest(times, time):
  idx = np.searchsorted(times, time, side='left')
  if idx > 0 and \
      (idx == len(times) or \
      math.fabs(time - times[idx-1]) < math.fabs(time - times[idx])):
    return idx - 1, times[idx - 1]
  else:
    return idx, times[idx]


def getInternalVariableSnapshots(uda_path, analytical_times, matID=0):

  # Get the internal variables
  eps_p_eq, epsdot_p_eq, backstress, porosity, damage, temperature, volume, times = getAllInternalVariables(
      uda_path, matID)

  # Create a clean list containing the data as functions of time
  idx_list = []
  time_list = []
  eps_p_list = []
  epsdot_p_list = []
  backstress_list = []
  porosity_list = []
  damage_list = []
  temp_list = []
  volume_list = []
  for ta in analytical_times:
    idx, time = find_nearest(times, ta)
    idx_list.append(idx)
    time_list.append(time)
    eps_p_list.append(eps_p_eq[idx])
    epsdot_p_list.append(epsdot_p_eq[idx])
    backstress_list.append(backstress[idx])
    porosity_list.append(porosity[idx])
    damage_list.append(damage[idx])
    temp_list.append(temperature[idx])
    volume_list.append(volume[idx])

  return idx_list, time_list, eps_p_list, epsdot_p_list, backstress_list, porosity_list, damage_list, temp_list, volume_list


#---------------------------------------------------------------------------------
# Read the internal state variables from the UDA file:
#---------------------------------------------------------------------------------
def getAllInternalVariables(uda_path, matID):

  # Get the internal variables
  times, eps_p_eq = get_pScalar(uda_path, "p.eqPlasticStrain", matID)
  times, epsdot_p_eq = get_pScalar(uda_path, "p.eqPlasticStrainRate", matID)
  times, backstress = get_pTensor(uda_path, "p.backStress", matID)

  # Get other variables that may be needed
  times, porosity = get_pScalar(uda_path, "p.porosity", matID)
  times, damage = get_pScalar(uda_path, "p.damage", matID)
  times, temperature = get_pScalar(uda_path, "p.temperature", matID)
  times, volume = get_pScalar(uda_path, "p.volume", matID)

  #print(times_list, eps_p_eq, epsdot_p_eq, backstress)

  return eps_p_eq, epsdot_p_eq, backstress, porosity, damage, temperature, volume, times


#---------------------------------------------------------------------------------
# Read a scalar variable from the UDA file:
#---------------------------------------------------------------------------------
def get_pScalar(uda_path, varname, matID=0):
  #Extract history
  print("Extracting", varname, "history...")
  args = [partextract_exe, "-partvar", varname, uda_path]
  F_scalar = tempfile.TemporaryFile()
  tmp = sub_proc.Popen(args, stdout=F_scalar, stderr=sub_proc.PIPE)
  dummy = tmp.wait()
  print('Done.')
  #Read file back in
  F_scalar.seek(0)
  times = []
  values = []
  FAIL_NAN = False
  for line in F_scalar:
    line = line.strip().split()
    times.append(float(line[0]))
    values.append(np.float64(line[4]))
    if np.isnan(values[-1]):
      FAIL_NAN = True
  if FAIL_NAN:
    print("ERROR: 'nan' encountered while retrieving", varname,
          ", will not plot correctly.")
  F_scalar.close()
  return times, values


#---------------------------------------------------------------------------------
# Read a 2-tensor variable from the UDA file:
#---------------------------------------------------------------------------------
def get_pTensor(uda_path, varname, matID=0):
  NAN_FAIL = False
  #Extract tensor history
  print("Extracting ", varname, " history...")
  args = [partextract_exe, "-partvar", varname, uda_path]
  print(args)
  F_tensor = tempfile.TemporaryFile()
  tmp = sub_proc.Popen(args, stdout=F_tensor, stderr=sub_proc.PIPE)
  dummy = tmp.wait()
  print('Done.')
  #Read file back in
  F_tensor.seek(0)
  times = []
  values = []
  for line in F_tensor:

    # If the first word in the line is Error then exit
    #print("Line = ", line)
    line = line.decode()
    first_word = line.partition(' ')[0]
    if first_word == "Error":
      print("**ERROR** partextract failed to read the", varname,
            "history data.")
      print("          It generated the following message:")
      print(line)
      sys.exit("Stopping program.")

    line = line.strip().split()
    times.append(float(line[0]))
    T11 = np.float64(line[4])
    T12 = np.float64(line[5])
    T13 = np.float64(line[6])
    T21 = np.float64(line[7])
    T22 = np.float64(line[8])
    T23 = np.float64(line[9])
    T31 = np.float64(line[10])
    T32 = np.float64(line[11])
    T33 = np.float64(line[12])
    values.append(np.array([[T11, T12, T13], [T21, T22, T23], [T31, T32, T33]]))
    for i in range(3):
      for j in range(3):
        if np.isnan(values[-1][i][j]):
          NAN_FAIL = True
  F_tensor.close()
  if NAN_FAIL:
    print("\nERROR: 'nan's found reading in", varname,
          ". Will not plot correctly")
  return times, values


def get_pStress(uda_path, matID=0):
  return get_pTensor(uda_path, "p.stress")


def get_pDeformationMeasure(uda_path, matID=0):
  return get_pTensor(uda_path, "p.deformationGradient")


def get_epsilons(uda_path):
  #Assumes no shear strains
  times, Fs = get_pDeformationMeasure(uda_path)
  epsils = []
  for F in Fs:
    epsils.append(
        np.array([[np.log(F[0][0]), 0, 0], [0, np.log(F[1][1]), 0],
                  [0, 0, np.log(F[2][2])]]))
  return times, epsils


def get_pEqStrain(uda_path, matID=0):
  return get_pScalar(uda_path, "p.eqStrain")


def get_pEqPlasticStrain(uda_path, matID=0):
  return get_pScalar(uda_path, "p.eqPlasticStrain")


def get_pEqPlasticStrainRate(uda_path, matID=0):
  return get_pScalar(uda_path, "p.eqPlasticStrainRate")


def get_totalEqStrain(uda_path):
  return get_pEqStrain(uda_path)


def get_defTable(uda_path, working_dir):
  #Determine the defTable file
  try:
    ups_file = os.path.abspath(uda_path) + '/input.xml.orig'
    F = open(ups_file, "r")
  except:
    ups_file = os.path.abspath(uda_path) + '/input.xml'
    F = open(ups_file, "r")

  def_file = ""
  for line in F:
    if '<prescribed_deformation_file>' in line and '</prescribed_deformation_file>' in line:
      def_file = line.split('<prescribed_deformation_file>')[1].split(
          '</prescribed_deformation_file>')[0].strip()
  F.close()

  #Assumes the input deck and uda share the same parent folder.
  def_file = os.path.join(working_dir, def_file)

  F = open(def_file, 'r')
  times = []
  Fs = []
  for line in F:
    line = line.strip().split()
    times.append(float(line[0]))
    Fs.append(
        np.array([[float(line[1]),
                   float(line[2]),
                   float(line[3])],
                  [float(line[4]),
                   float(line[5]),
                   float(line[6])],
                  [float(line[7]),
                   float(line[8]),
                   float(line[9])]]))
  F.close()
  return times, Fs


def exp_fmt(x, loc):
  tmp = format(x, '1.2e').split('e')
  lead = tmp[0]
  exp = str(int(tmp[1]))
  if exp == '0' and lead == '0.00':
    return r'$\mathbf{0.00}$'
  else:
    if int(exp) < 10 and int(exp) > 0:
      exp = '+0' + exp
    elif int(exp) > -10 and int(exp) < 0:
      exp = '-0' + exp.split('-')[1]
    elif int(exp) > 10:
      exp = '+' + exp
    return r'$\mathbf{' + lead + r'\cdot{}10^{' + exp + '}}$'


#------------------------------------------------------------------------
# Plot x vs y
#------------------------------------------------------------------------
def plotEqShearMeanStress(xs,
                          ys,
                          idx_list,
                          color_list,
                          eps_p_eq_list,
                          compression='negative',
                          Xlims=False,
                          Ylims=False,
                          GRID=True):

  ax1 = plt.subplot(111)
  if (compression == 'positive'):
    xs = list(map(lambda x: -x, xs))
    ys = list(map(lambda y: -y, ys))

  x_data = np.array(xs)
  y_data = np.array(ys)

  max_idx = max(idx_list)
  for jj, color in enumerate(color_list):
    eps_p_eq_str = str(eps_p_eq_list[jj])
    label_str = 'Eq. plastic strain = ' + eps_p_eq_str
    start_idx = idx_list[jj]
    if (start_idx < max_idx):
      end_idx = idx_list[jj + 1]
      x_arrow = x_data[start_idx:end_idx - 1]
      y_arrow = y_data[start_idx:end_idx - 1]
      u_arrow = x_data[start_idx + 1:end_idx] - x_arrow
      v_arrow = y_data[start_idx + 1:end_idx] - y_arrow
      plt.quiver(x_arrow,
                 y_arrow,
                 u_arrow,
                 v_arrow,
                 scale_units='xy',
                 angles='xy',
                 scale=1.5,
                 width=0.001,
                 headwidth=7,
                 headlength=10,
                 linewidth=0.5,
                 edgecolor='b',
                 color='k')
      plt.plot(x_data[start_idx:end_idx], y_data[start_idx:end_idx], '-', \
        color=color,label=label_str)

  plt.xlabel(str_to_mathbf('Mean Stress, p (Pa)'))
  plt.ylabel(str_to_mathbf('Equivalent Shear Stress, q, (Pa)'))

  formatter_int = ticker.FormatStrFormatter('$\mathbf{%g}$')
  formatter_exp = ticker.FuncFormatter(exp_fmt)

  #ax1.xaxis.set_major_formatter(formatter_exp)
  #ax1.yaxis.set_major_formatter(formatter_exp)
  ax1.xaxis.set_major_formatter(formatter_int)
  ax1.yaxis.set_major_formatter(formatter_int)

  if Xlims:
    ax1.set_xlim(Xlims[0], Xlims[1])
  if Ylims:
    ax1.set_ylim(Ylims[0], Ylims[1])
  if GRID:
    plt.grid(True)

  return ax1


def eqShear_vs_meanStress(xs,
                          ys,
                          compression='negative',
                          Xlims=False,
                          Ylims=False,
                          LINE_LABEL='Vaango',
                          GRID=True):

  ax1 = plt.subplot(111)
  if (compression == 'negative'):
    x_data = np.array(xs)
    y_data = np.array(ys)
    plt.quiver(x_data[:-1],
               y_data[:-1],
               x_data[1:] - x_data[:-1],
               y_data[1:] - y_data[:-1],
               scale_units='xy',
               angles='xy',
               scale=1.0,
               width=0.001,
               headwidth=7,
               headlength=10,
               linewidth=0.5,
               edgecolor='b',
               color='k')
    plt.plot(x_data, y_data, '-r', label=LINE_LABEL)
  else:
    xs = list(map(lambda x: -x, xs))
    ys = list(map(lambda y: -y, ys))
    x_data = np.array(xs)
    y_data = np.array(ys)
    plt.quiver(x_data[:-1],
               y_data[:-1],
               x_data[1:] - x_data[:-1],
               y_data[1:] - y_data[:-1],
               scale_units='xy',
               angles='xy',
               scale=1.0,
               width=0.001,
               headwidth=7,
               headlength=10,
               linewidth=0.5,
               edgecolor='b',
               color='k')
    plt.plot(x_data, y_data, '-r', label=LINE_LABEL)

  plt.xlabel(str_to_mathbf('Mean Stress, p (Pa)'))
  plt.ylabel(str_to_mathbf('Equivalent Shear Stress, q, (Pa)'))

  formatter_int = ticker.FormatStrFormatter('$\mathbf{%g}$')
  formatter_exp = ticker.FuncFormatter(exp_fmt)

  #ax1.xaxis.set_major_formatter(formatter_exp)
  #ax1.yaxis.set_major_formatter(formatter_exp)
  ax1.xaxis.set_major_formatter(formatter_int)
  ax1.yaxis.set_major_formatter(formatter_int)

  if Xlims:
    ax1.set_xlim(Xlims[0], Xlims[1])
  if Ylims:
    ax1.set_ylim(Ylims[0], Ylims[1])
  if GRID:
    plt.grid(True)

  return ax1


def I1_to_zbar(I1s):
  sqrt_3 = np.sqrt(3.0)
  if type(I1s) in [list, np.ndarray]:
    zbars = []
    for I1 in I1s:
      zbars.append(-I1 / sqrt_3)
    return zbars
  elif type(I1s) in [int, float, np.float64]:
    return -I1s / sqrt_3
  else:
    print(
        '\nERROR: cannot compute zbar from I1. Invalid type.\n\ttype(I1)\t:\t',
        type(I1s))
    return None


#-----------------------------------------------------------------------------------
# Read the simulation stress data and compute p,q (converts to Pa)
#-----------------------------------------------------------------------------------
def readSimStressData(uda_path, matID=0):

  NAN_FAIL = False

  #Extract stress history
  time_sim, sigma_sim = get_pTensor(uda_path, "p.stress")
  #print("stress = ", sigma_sim)

  sigma_a_sim = [sigma[0][0] for sigma in sigma_sim]
  sigma_r_sim = [sigma[1][1] for sigma in sigma_sim]
  sigma_ar_sim = [sigma[0][1] for sigma in sigma_sim]

  # Compute I1 and J2
  I1_sim = [sigma_I1(sigma) for sigma in sigma_sim]
  J2_sim = [sigma_J2(sigma) for sigma in sigma_sim]
  J3_sim = [sigma_J3(sigma) for sigma in sigma_sim]

  # Compute p, q
  pp_sim = list(map(lambda I1: I1 / 3, I1_sim))
  qq_sim = list(
      map(lambda J2, J3: np.sign(J3) * np.sqrt(3 * J2), J2_sim, J3_sim))

  return time_sim, sigma_sim, sigma_a_sim, sigma_r_sim, sigma_ar_sim, pp_sim, qq_sim


#---------------------------------------------------------------------------------
# Get the closest data at a set of given time points
#---------------------------------------------------------------------------------
def getDataTimeSnapshots(time_snapshots, times, data):

  # Create a clean list containing the data as functions of time
  time_list = []
  val_list = []
  for ta in time_snapshots:
    idx, time = find_nearest(times, ta)
    time_list.append(time)
    val_list.append(data[idx])

  return time_list, val_list


#-----------------------------------------------------------------------------------
# Plot the sim data as a function of stress vs strain
#-----------------------------------------------------------------------------------
def plotSimDataSigmaEps(fig,
                        time_snapshots,
                        time_sim,
                        pp_sim,
                        ev_e_sim,
                        ev_p_sim,
                        compression='negative'):

  # Get snapshots from data
  time_snap, pp_snap = getDataTimeSnapshots(time_snapshots, time_sim, pp_sim)

  # Activate the figure
  plt.figure(fig.number)

  ev_sim = list(map(lambda ev_e, ev_p: ev_e + ev_p, ev_e_sim, ev_p_sim))

  # Plot sigma_a vs. time
  if (compression == 'positive'):
    ev_data = list(map(lambda p: -p, ev_sim))
    pp_data = list(map(lambda p: -p, pp_sim))
    plt.plot(ev_data, pp_data, '--r', label='Simulation')
  else:
    plt.plot(ev_sim, pp_sim, '--r', label='Simulation')

  return time_snap, pp_snap


#-----------------------------------------------------------------------------------
# Plot the sim data as a function of time (assume Pa)
#-----------------------------------------------------------------------------------
def plotSimDataSigmaTime(fig,
                         time_snapshots,
                         time_sim,
                         sigma_a_sim,
                         sigma_r_sim,
                         sigma_ar_sim,
                         labelxx='$\sigma_a$ (sim)',
                         labelyy='$\sigma_r$ (sim)',
                         labelxy='$\sigma_{ar}$ (sim)',
                         compression='negative'):

  # Get snapshots from data
  time_snap, sigma_a_snap = getDataTimeSnapshots(time_snapshots, time_sim,
                                                 sigma_a_sim)
  time_snap, sigma_r_snap = getDataTimeSnapshots(time_snapshots, time_sim,
                                                 sigma_r_sim)
  time_snap, sigma_ar_snap = getDataTimeSnapshots(time_snapshots, time_sim,
                                                  sigma_ar_sim)

  # Activate the figure
  plt.figure(fig.number)

  # Plot sigma_a vs. time
  if (compression == 'positive'):
    sigma_a_data = list(map(lambda p: -p, sigma_a_sim))
    sigma_r_data = list(map(lambda p: -p, sigma_r_sim))
    sigma_ar_data = list(map(lambda p: -p, sigma_ar_sim))
    plt.plot(time_sim, sigma_a_data, '--r', label=labelxx)
    plt.plot(time_sim, sigma_r_data, '--b', label=labelyy)
    plt.plot(time_sim, sigma_ar_data, '--g', label=labelxy)
  else:
    plt.plot(time_sim, sigma_a_sim, '--r', label=labelxx)
    plt.plot(time_sim, sigma_r_sim, '--b', label=labelyy)
    plt.plot(time_sim, sigma_ar_sim, '--g', label=labelxy)

  # Plot filled circles at time snapshots
  for ii in range(0, len(time_snap)):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(time_snap))
    if (compression == 'positive'):
      plt.plot(time_snap[ii], -sigma_a_snap[ii], 'o', color=plt_color)
      plt.plot(time_snap[ii], -sigma_r_snap[ii], 'o', color=plt_color)
    else:
      plt.plot(time_snap[ii], sigma_a_snap[ii], 'o', color=plt_color)
      plt.plot(time_snap[ii], sigma_r_snap[ii], 'o', color=plt_color)
      #plt.plot(time_snap[ii], sigma_ar_snap[ii], 'o', color=plt_color)

  return time_snap, sigma_a_snap, sigma_r_snap, sigma_ar_snap


#-----------------------------------------------------------------------------------
# Plot the sim data (pq) as a function of time (aasume Pa)
#-----------------------------------------------------------------------------------
def plotSimDataPQTime(fig,
                      time_snapshots,
                      time_sim,
                      p_sim,
                      q_sim,
                      compression='negative'):

  # Get snapshots from sim data
  time_snap, p_snap = getDataTimeSnapshots(time_snapshots, time_sim, p_sim)
  time_snap, q_snap = getDataTimeSnapshots(time_snapshots, time_sim, q_sim)

  # Activate the figure
  plt.figure(fig.number)

  # Plot sigma_a vs. time
  if (compression == 'positive'):
    p_data = list(map(lambda p: -p, p_sim))
    q_data = list(map(lambda q: -q, q_sim))
    plt.plot(time_sim, p_data, '--r', label='$p$ (sim)')
    plt.plot(time_sim, q_data, '--b', label='$q$ (sim)')
  else:
    plt.plot(time_sim, p_sim, '--r', label='$p$ (sim)')
    plt.plot(time_sim, q_sim, '--b', label='$q$ (sim)')

  # Plot filled circles at time snapshots
  for ii in range(0, len(time_snap)):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(time_snap))
    if (compression == 'positive'):
      plt.plot(time_snap[ii], -p_snap[ii], 'o', color=plt_color)
      plt.plot(time_snap[ii], -q_snap[ii], 'o', color=plt_color)
    else:
      plt.plot(time_snap[ii], p_snap[ii], 'o', color=plt_color)
      plt.plot(time_snap[ii], q_snap[ii], 'o', color=plt_color)

  return time_snap, p_snap, q_snap
