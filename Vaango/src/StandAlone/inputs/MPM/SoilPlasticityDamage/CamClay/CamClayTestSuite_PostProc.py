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
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker

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
  plt.gcf().set_size_inches(size[0], size[1])
  #save at speciified resolution
  plt.savefig(name + '.png', bbox_inches=0, dpi=plt.rcParams['figure.dpi'])


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


def str_to_mathbf(string):
  #Only works with single spaces no leading space
  string = string.split()
  return_string = ''
  for elem in string:
    elem = r'$\mathbf{' + elem + '}$'
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


def get_ps_and_qs(sigmas):
  ps = []
  qs = []
  for sigma in sigmas:
    qs.append(sign(sqrtThree * np.sqrt(sigma_J2(sigma)), sigma_J3(sigma)))
    ps.append(sigma_I1(sigma) / 3.0)
  return ps, qs

def get_pScalar(uda_path, varname):
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
    print(
        "ERROR: 'nan' encountered while retrieving", varname, ", will not plot correctly."
    )
  F_scalar.close()
  return times, values

def get_pTensor(uda_path, varname):
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
      print("**ERROR** partextract failed to read the", varname, "history data.")
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
    print("\nERROR: 'nan's found reading in", varname, ". Will not plot correctly")
  return times, values

def get_epsilons(uda_path):
  #Assumes no shear strains
  times, Fs = get_pTensor(uda_path, "p.deformationGradient")
  epsils = []
  for F in Fs:
    epsils.append(
        np.array([[np.log(F[0][0]), 0, 0], [0, np.log(F[1][1]), 0],
                  [0, 0, np.log(F[2][2])]]))
  return times, epsils



def get_defTable(uda_path, working_dir):
  #Determine the defTable file
  try:
    ups_file = os.path.abspath(uda_path) + '/input.xml.orig'
    F = open(ups_file, "r")
  except:
    ups_file = os.path.abspath(uda_path) + '/input.xml'
    F = open(ups_file, "r")
  for line in F:
    if '<prescribed_deformation_file>' in line and '</prescribed_deformation_file>' in line:
      def_file = line.split('<prescribed_deformation_file>')[1].split(
          '</prescribed_deformation_file>')[0].strip()
  F.close()
  #Assumes the input deck and uda share the same parent folder.
  def_file = working_dir + '/' + def_file
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


def eqShear_vs_meanStress(xs,
                          ys,
                          Xlims=False,
                          Ylims=False,
                          LINE_LABEL='Vaango',
                          GRID=True):

  ax1 = plt.subplot(111)
  plt.plot(np.array(xs), np.array(ys), '-r', label=LINE_LABEL)

  plt.xlabel(str_to_mathbf('Mean Stress, p (MPa)'))
  plt.ylabel(str_to_mathbf('Equivalent Shear Stress, q, (MPa)'))

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


def get_yield_surface(uda_path):
  try:
    ups_file = os.path.abspath(uda_path) + '/input.xml.orig'
    F_ups = open(ups_file, "r")
  except:
    ups_file = os.path.abspath(uda_path) + '/input.xml'
    F_ups = open(ups_file, "r")

  check_lines = False
  already_read = False
  material_dict = {}
  for line in F_ups:
    if '<constitutive_model' in line and 'type' in line and '"camclay"' in line and not (
        already_read):
      check_lines = True
    if check_lines and not (already_read):
      if '<p0>' in line:
        material_dict['p0'] = float(
            line.split('<p0>')[1].split('</p0>')[0].strip())
      if '<alpha>' in line:
        material_dict['alpha'] = float(
            line.split('<alpha>')[1].split('</alpha>')[0].strip())
      if '<kappatilde>' in line:
        material_dict['kappatilde'] = float(
            line.split('<kappatilde>')[1].split('</kappatilde>')[0].strip())
      if '<epse_v0>' in line:
        material_dict['epse_v0'] = float(
            line.split('<epse_v0>')[1].split('</epse_v0>')[0].strip())
      if '<mu0>' in line:
        material_dict['mu0'] = float(
            line.split('<mu0>')[1].split('</mu0>')[0].strip())
      if '<M>' in line:
        material_dict['M'] = float(
            line.split('<M>')[1].split('</M>')[0].strip())
      if '<pc0>' in line:
        material_dict['pc0'] = float(
            line.split('<pc0>')[1].split('</pc0>')[0].strip())
      if '<lambdatilde>' in line:
        material_dict['lambdatilde'] = float(
            line.split('<lambdatilde>')[1].split('</lambdatilde>')[0].strip())
      if '</constitutive_model>' in line:
        already_read = True
        check_lines = False
  F_ups.close()
  PRINTOUT = False
  if PRINTOUT:
    print('--Material Specification--')
    for key in material_dict:
      print(key, ':', material_dict[key])

  #tmp_string = r'$\mathbf{\underline{Material}}$'+' '+r'$\mathbf{\underline{Properties:}}$'+'\n'
  tmp_string = r'$\mathbf{\underline{Material\phantom{1}Properties:}}$' + '\n'

  sorted_list_tuples = sorted(material_dict.items())
  #print(sorted_list_tuples)
  sorted_dict = {}
  for key, value in sorted_list_tuples:
    if '_' in key:
      tmp = key.split('_')
      tmp = str_to_mathbf(tmp[0] + '_' + '{' + tmp[1] + '}')
      tmp_string += tmp + str_to_mathbf(' = ') + str_to_mathbf(
          format(value, '1.3e')) + '\n'
    else:
      tmp_string += str_to_mathbf(key + ' = ' + format(value, '1.3e')) + '\n'
  material_dict['material string'] = tmp_string[0:-1]
  if PRINTOUT:
    print(tmp_string)
  return material_dict


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


def J2VM(epsil_dot, dt, sig_Beg, K, G, tau_y):
  #J2 plasticity Von misses material model for 3D
  #Inputs: epsil_dot, dt, sig_Beg, K, G, tau_y
  #Outputs: epsil_Elastic_dot, epsil_Plastic_dot, sig_End
  #Initialize the trial stress state
  sig_Trial = sig_Beg + (
      (2 * G * sigma_dev(epsil_dot)) + 3 * K * sigma_iso(epsil_dot)) * dt
  #Determine if this is below, on, or above the yeild surface
  test = sigma_mag(sigma_dev(sig_Trial)) / (np.sqrt(2.0) * tau_y)
  if test <= 1:
    #Stress state is elastic
    sig_End = sig_Trial
    epsil_Plastic_dot = np.zeros((3, 3))
    epsil_Elastic_dot = epsil_dot
  elif test > 1:
    #Stress state elastic-plastic
    sig_End = (sigma_dev(sig_Trial) / test)  #+sigma_iso(sig_Trial)
    #Evaluate the consistent stress rate
    #sig_dot = (sig_End-sig_Beg)/test
    #Apply hookes law to get the elastic strain rate
    #epsil_Elastic_dot = sigma_dev(sig_dot)/(2*G)# + sigma_iso(sig_dot)/(3*K)
    #Apply strain rate decomposition relationship to get plastic strain rate
    #epsil_Plastic_dot = epsil_dot-epsil_Elastic_dot
  #Determine the equivalent stress and equivalent plastic strain rate
  #sig_Eq = np.sqrt(3/2)*sigma_mag(sigma_dev(sig_End))
  #epsil_Plastic_dot_Eq = np.sqrt(3/2)*sigma_mag(sigma_dev(epsil_Plastic_dot))
  #ans={'Elastic dot':epsil_Elastic_dot,'Plastic dot':epsil_Plastic_dot,'Stress State':sig_End}
  return sig_End


def defTable_to_J2Solution(def_times,
                           Fs,
                           bulk_mod,
                           shear_mod,
                           tau_yield,
                           num_substeps=1000):
  #Assumes:
  print('Solving for analytical solution...')
  analytical_epsils = [np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]
  analytical_sigmas = [np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]
  analytical_times = [def_times[0]]

  epsils = []
  for F in Fs:
    epsils.append(
        np.array([[np.log(sum(F[0])), 0, 0], [0, np.log(sum(F[1])), 0],
                  [0, 0, np.log(sum(F[2]))]]))

  for leg in range(len(def_times) - 1):
    t_start = def_times[leg]
    leg_delT = def_times[leg + 1] - t_start
    leg_sub_delT = float(leg_delT) / float(num_substeps)
    leg_del_epsil = (epsils[leg + 1] - epsils[leg])
    leg_epsil_dot = leg_del_epsil / leg_delT
    for i in range(num_substeps):
      t_now = t_start + float(i) * leg_sub_delT
      analytical_times.append(t_now)
      analytical_sigmas.append(
          J2VM(leg_epsil_dot, leg_sub_delT, analytical_sigmas[-1], bulk_mod,
               shear_mod, tau_yield))
      analytical_epsils.append(analytical_epsils[-1] +
                               (leg_epsil_dot * leg_sub_delT))
  analytical_epsils.append(analytical_epsils[-1] +
                           (leg_epsil_dot * leg_sub_delT))
  analytical_sigmas.append(
      J2VM(leg_epsil_dot, leg_sub_delT, analytical_sigmas[-1], bulk_mod,
           shear_mod, tau_yield))
  analytical_times.append(def_times[-1])
  print('Done.')
  return analytical_times, analytical_sigmas, analytical_epsils


def plot_yield_surface(uda_path, PLOT_TYPE='J2_vs_I1'):
  num_points = 500
  material_dict = get_yield_surface(uda_path)
  p0 = material_dict['p0']
  alpha = material_dict['alpha']
  kappatilde = material_dict['kappatilde']
  epse_v0 = material_dict['epse_v0']
  mu0 = material_dict['mu0']
  M = material_dict['M']
  pc0 = material_dict['pc0']
  lambdatilde = material_dict['lambdatilde']

  # Compute p/q
  ps = np.linspace(p0, pc0, num_points)
  qs = [M * np.sqrt(np.abs(p * (pc0 - p))) for p in ps]

  # Convert to I1/J2
  I1s = [3.0 * p for p in ps]
  J2s = [q * q / 3 for q in qs]

  # Convert to r/z
  rs = [I1 / np.sqrt(3) for I1 in I1s]
  zs = [np.sqrt(2.0 * J2) for J2 in J2s]

  PLOT = True
  xs = []
  ys = []
  if PLOT_TYPE == 'J2_vs_I1':
    xs = np.array(I1s) * 1.0e-6
    ys = np.array(J2s) * 1.0e-12
  elif PLOT_TYPE == 'sqrtJ2_vs_I1':
    xs = np.array(I1s) * 1.0e-6
    ys = np.sqrt(np.array(J2s)) * 1.0e-6
  elif PLOT_TYPE == 'r_vs_z':
    xs = np.array(rs) * 1.0e-6
    ys = np.array(zs) * 1.0e-6
  elif PLOT_TYPE == 'q_vs_I1':
    xs = np.array(I1s) * 1.0e-6
    ys = np.array(qs) * 1.0e-6
  elif PLOT_TYPE == 'q_vs_p':
    xs = np.array(ps) * 1.0e-6
    ys = np.array(qs) * 1.0e-6
  else:
    PLOT = False
    print(
        '\nError: invalid plot type specified for initial yield surface plot.\n\tPLOT_TYPE:',
        PLOT_TYPE)
  if PLOT:
    plt.plot(xs,
             ys,
             '--k',
             linewidth=lineWidth + 1,
             label='Initial Yield Surface')
    plt.plot(xs, -ys, '--k', linewidth=lineWidth + 1)


def test_yield_surface(uda_path):
  plot_yield_surface(uda_path, 'J2_vs_I1')
  plt.show()
  plot_yield_surface(uda_path, 'sqrtJ2_vs_I1')
  plt.show()
  plot_yield_surface(uda_path, 'r_vs_z')
  plt.show()
  plot_yield_surface(uda_path, 'q_vs_I1')
  plt.show()
  plot_yield_surface(uda_path, 'q_vs_p')
  plt.show()


### ----------
#  Test Methods Below
### ----------
def test01_postProc(uda_path, save_path, **kwargs):
  print("CamClay Test: 01 - Uniaxial Compression With Rotation")

  # Read the simulation data
  times, sigmas, sigma_a_sim, sigma_r_sim, sigma_ar_sim, pp_sim, qq_sim = readSimStressData(
      uda_path)

  #print(times)
  #print(sigmas)

  # Get the model parameters
  material_dict = get_yield_surface(uda_path)
  param_text = material_dict['material string']

  # Set up time points
  #analytical_times = [0.0, 0.1, 0.2, 0.5, 0.7, 1.0]
  analytical_times = np.linspace(0.0, 0.001, 15)

  # Get snapshots of pq data (sim)
  t_sim_snap, p_sim_snap = getDataTimeSnapshots(analytical_times, times, pp_sim)
  t_sim_snap, q_sim_snap = getDataTimeSnapshots(analytical_times, times, qq_sim)

  # Find the plot limits
  Sxx = []
  Syy = []
  for sigma in sigmas:
    Sxx.append(sigma[0][0])
    Syy.append(sigma[1][1])

  # Find min/max values
  Sxx_min = min(Sxx)
  Syy_min = min(Syy)
  Sxx_max = max(Sxx)
  Syy_max = max(Syy)
  print("Sxx_min = ", Sxx_min)
  print("Sxx_max = ", Sxx_max)
  print("Syy_min = ", Syy_min)
  print("Syy_max = ", Syy_max)
  Sxx_tick_int = (Sxx_max - Sxx_min) / 4
  Syy_tick_int = (Syy_max - Syy_min) / 4
  #Sxx_min_f = math.floor(Sxx_min)
  #Syy_min_f = math.floor(Syy_min)
  #Sxx_max_c = math.ceil(Sxx_max)
  #Syy_max_c = math.ceil(Syy_max)
  Sxx_ticks = [
      round(Sxx_max + Sxx_tick_int, 1),
      round(Sxx_max, 1),
      round(Sxx_max - Sxx_tick_int, 1),
      round(Sxx_max - 2 * Sxx_tick_int, 1),
      round(Sxx_max - 3 * Sxx_tick_int, 1),
      round(Sxx_max - 4 * Sxx_tick_int, 1)
  ]
  Syy_ticks = [
      round(Syy_max + Syy_tick_int, 1),
      round(Syy_max, 1),
      round(Syy_max - Syy_tick_int, 1),
      round(Syy_max - 2 * Syy_tick_int, 1),
      round(Syy_max - 3 * Syy_tick_int, 1),
      round(Syy_max - 4 * Syy_tick_int, 1)
  ]

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  plt.figure(1)
  plt.clf()
  plt.subplots_adjust(right=0.75)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  #Syy
  ax2 = plt.subplot(212)
  #without rotation
  #plt.plot([0,1],[0,0],'-b')
  #simulation results
  plt.plot(times, Syy, '-r')
  #guide line
  #plt.plot([0,1],[Syy_max,Syy_min],'--g')
  #labels and limits
  #ax2.set_xlim(0, 1);
  ax2.set_ylim(Syy_min, Syy_max + Syy_tick_int)
  ax2.set_yticks(Syy_ticks)
  plt.grid(True)
  ax2.xaxis.set_major_formatter(formatter)
  ax2.yaxis.set_major_formatter(formatter)
  plt.ylabel(str_to_mathbf('\sigma_{yy} (MPa)'))
  plt.xlabel(str_to_mathbf('Time (s)'))
  #Sxx
  ax1 = plt.subplot(211, sharex=ax2, sharey=ax2)
  plt.setp(ax1.get_xticklabels(), visible=False)
  #without rotation
  #plt.plot([0,1],[Sxx_max, Sxx_min],'-b',label='No rotation')
  #simulation results
  plt.plot(times, Sxx, '-r', label='Vaango')
  #guide lines
  #plt.plot([0,1],[0,0],'--g',label='Guide lines')
  #labels
  ax1.set_ylim(Sxx_min, Sxx_max + Sxx_tick_int)
  #ax1.set_xlim(0, 1)
  plt.grid(True)
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  ax1.set_yticks(Sxx_ticks)
  plt.ylabel(str_to_mathbf('\sigma_{xx} (MPa)'))
  plt.title('CamClay 01:\nUniaxial Strain Compression')
  plt.legend()
  savePNG(save_path + '/Test01_verificationPlot', '1280x960')
  if SHOW_ON_MAKE:
    plt.show()

  #----------------------------------------------------------------
  # Plot the yield surface for test1
  #----------------------------------------------------------------

  # Set up figure
  fig2 = plt.figure(2)
  plt.clf()
  plt.subplots_adjust(right=0.75)
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')

  # Change sign
  pp_sim = [-p for p in pp_sim]
  qq_sim = [-q for q in qq_sim]
  p_sim_snap = [-p for p in p_sim_snap]
  q_sim_snap = [-q for q in q_sim_snap]

  # Plot p vs. q simulation results
  eqShear_vs_meanStress(pp_sim, qq_sim)

  # Plot filled circles at time snapshots
  for ii in range(0, len(t_sim_snap) - 1):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(t_sim_snap))
    plt.plot(p_sim_snap[ii], q_sim_snap[ii], 'o', color=plt_color)

  # Plot yield surfaces
  pMin, qMax = plotPQYieldSurfaceSim(uda_path,
                                     analytical_times)

  plt.title(
      'CamClay: Uniaxial strain compression\n Yield surface evolution'
  )
  savePNG(save_path + '/Test01_yield_surface', '1280x960')
  #plt.show()

  #---------------------------------------------------------------------------------
  # Plot simulation data as a function of time
  fig3 = plt.figure(3)
  plt.clf()
  plt.subplots_adjust(right=0.75)
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  plotSimDataSigmaTime(fig3, analytical_times, times, sigma_a_sim, sigma_r_sim,
                       sigma_ar_sim)
  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  plt.xlabel(str_to_mathbf('Time (sec)'))
  plt.ylabel(str_to_mathbf('Stress (MPa)'))
  plt.grid(True)
  plt.legend(loc='best', prop={'size': 10})
  savePNG(save_path + '/Test01_sigma_time', '1280x960')

  fig4 = plt.figure(4)
  plt.clf()
  plt.subplots_adjust(right=0.75)
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  plotSimDataPQTime(fig4, analytical_times, times, pp_sim, qq_sim)
  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  plt.xlabel(str_to_mathbf('Time (sec)'))
  plt.ylabel(str_to_mathbf('Stress invariants (p/q) (MPa)'))
  plt.grid(True)
  plt.legend(loc='best', prop={'size': 10})
  savePNG(save_path + '/Test01_pq_time', '1280x960')

  plt.show()


def test02_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 02 - Vertex Treatment")
  times, sigmas = get_pTensor(uda_path, "p.stress")
  time_list, pcs = get_pScalar(uda_path, "p.p_c")
  ps_unscaled, qs_unscaled = get_ps_and_qs(sigmas)
  ps = []
  qs = []
  for val in ps_unscaled:
    ps.append(val * 1.0e-6)

  for val in qs_unscaled:
    qs.append(val * 1.0e-6)

  Sxx = []
  Syy = []
  Szz = []
  for sigma in sigmas:
    Sxx.append(sigma[0][0] * 1.0e-6)
    Syy.append(sigma[1][1] * 1.0e-6)
    Szz.append(sigma[2][2] * 1.0e-6)

  # Find min/max values
  ps_min = min(ps)
  ps_max = max(ps)
  qs_min = min(qs)
  qs_max = max(qs)
  Sxx_min = min(Sxx)
  Syy_min = min(Syy)
  Szz_min = min(Szz)
  Sxx_max = max(Sxx)
  Syy_max = max(Syy)
  Szz_max = max(Szz)
  print("Sxx_min = ", Sxx_min)
  print("Sxx_max = ", Sxx_max)
  print("Syy_min = ", Syy_min)
  print("Syy_max = ", Syy_max)
  print("Szz_min = ", Szz_min)
  print("Szz_max = ", Szz_max)
  Sxx_tick_int = (Sxx_max - Sxx_min) / 4
  Syy_tick_int = (Syy_max - Syy_min) / 4
  Szz_tick_int = (Szz_max - Szz_min) / 4
  Sxx_ticks = [
      round(Sxx_max + Sxx_tick_int, 1),
      round(Sxx_max, 1),
      round(Sxx_max - Sxx_tick_int, 1),
      round(Sxx_max - 2 * Sxx_tick_int, 1),
      round(Sxx_max - 3 * Sxx_tick_int, 1),
      round(Sxx_max - 4 * Sxx_tick_int, 1)
  ]
  Syy_ticks = [
      round(Syy_max + Syy_tick_int, 1),
      round(Syy_max, 1),
      round(Syy_max - Syy_tick_int, 1),
      round(Syy_max - 2 * Syy_tick_int, 1),
      round(Syy_max - 3 * Syy_tick_int, 1),
      round(Syy_max - 4 * Syy_tick_int, 1)
  ]
  Szz_ticks = [
      round(Szz_max + Szz_tick_int, 1),
      round(Szz_max, 1),
      round(Szz_max - Szz_tick_int, 1),
      round(Szz_max - 2 * Szz_tick_int, 1),
      round(Szz_max - 3 * Szz_tick_int, 1),
      round(Szz_max - 4 * Szz_tick_int, 1)
  ]

  analytical_times = [0,1,3.0/2.0,2.0,5.0/2.0,3.0]
  pc_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        pc_list.append(pcs[ii])

  #print("pc = ", pc_list)
  pc_list_new = list(sorted(set(pc_list)))

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  ##Plot a

  plt.figure(1)
  plt.clf()
  plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  material_dict = get_yield_surface(uda_path)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')

  ps = [-p for p in ps]
  qs = [-q for q in qs]
  eqShear_vs_meanStress(ps, qs)

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('CamClay 02:\nVertex Treatment (plot a)')
  #plt.legend()
  savePNG(save_path + '/Test02_verificationPlot_a', '1280x960')

  ##Plot b
  plt.figure(2)
  plt.clf()
  plt.subplots_adjust(right=0.75)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  endT = max(times)
  #Sigma zz
  ax3 = plt.subplot(313)
  plt.plot(times, np.array(Szz), '-r')
  #Add Yield Surface
  #Add Analytical
  plt.xlabel(str_to_mathbf('Time (s)'))
  plt.ylabel(str_to_mathbf('\sigma_{zz} (MPa)'))
  ax3.yaxis.set_major_formatter(formatter)
  ax3.set_xlim(0, endT)
  #ax3.set_ylim(-300,100)
  #ax3.set_yticks([-300,-200,-100,0,100])
  plt.grid(True)
  #Sigma xx
  ax1 = plt.subplot(311, sharex=ax3)
  plt.plot(times, np.array(Sxx), '-r', label='Vaango')
  #Add Yield Surface
  #Add Analytical
  plt.legend()
  plt.setp(ax1.get_xticklabels(), visible=False)
  plt.ylabel(str_to_mathbf('\sigma_{xx} (MPa)'))
  plt.title('CamClay 02:\nVertex Treatment (plot b)')
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  ax1.set_xlim(0, endT)
  #ax1.set_ylim(-400,100)
  #ax1.set_yticks([-400,-300,-200,-100,0,100])
  plt.grid(True)
  #Sigma yy
  ax2 = plt.subplot(312, sharex=ax3)
  plt.plot(times, np.array(Syy), '-r')
  #Add Yield Surface
  #Add Analytical
  plt.setp(ax2.get_xticklabels(), visible=False)
  plt.ylabel(str_to_mathbf('\sigma_{yy} (MPa)'))
  ax2.yaxis.set_major_formatter(formatter)
  ax2.set_xlim(0, endT)
  #ax2.set_ylim(-300,100)
  #ax2.set_yticks([-300,-200,-100,0,100])
  plt.grid(True)
  savePNG(save_path + '/Test02_verificationPlot_b', '1280x960')

  plt.show()


def test03_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 03 - Uniaxial Strain Without Hardening")
  times, sigmas = get_pTensor(uda_path, "p.stress")
  times, pcs = get_pScalar(uda_path, "p.p_c")
  ps_unscaled, qs_unscaled = get_ps_and_qs(sigmas)
  material_dict = get_yield_surface(uda_path)

  # Scale the data
  ps = [-p*1.0e-6 for p in ps_unscaled]
  qs = [-q*1.0e-6 for q in qs_unscaled]

  # Find min/max values
  ps_min = min(ps)
  ps_max = max(ps)
  qs_min = min(qs)
  qs_max = max(qs)
  print("ps_min = ", ps_min)
  print("ps_max = ", ps_max)
  print("qs_min = ", qs_min)
  print("qs_max = ", qs_max)

  analytical_times = [0.0, 1.0, 2.0]

  pc_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        pc_list.append(pcs[ii])

  #print("pc = ", pc_list)
  pc_list_new = list(sorted(set(pc_list)))

  ###PLOTTING
  Xlims = (ps_min, ps_max)
  Ylims = (qs_min, -qs_min)
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  material_dict = get_yield_surface(uda_path)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  #eqShear_vs_meanStress(ps,qs,Xlims,Ylims,)
  eqShear_vs_meanStress(ps, qs)
  plt.title('CamClay 03:\nUniaxial Strain Without Hardening')

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  #plt.legend()
  savePNG(save_path + '/Test03_verificationPlot', '1280x960')
  plt.show()


def test04_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 04 - Curved Yield Surface")
  times, sigmas = get_pTensor(uda_path, "p.stress")
  times, pcs = get_pScalar(uda_path, "p.p_c")
  ps_unscaled, qs_unscaled = get_ps_and_qs(sigmas)

  # Scale the data
  ps = [-p*1.0e-6 for p in ps_unscaled]
  qs = [-q*1.0e-6 for q in qs_unscaled]

  # Find min/max values
  ps_min = min(ps)
  ps_max = max(ps)
  qs_min = min(qs)
  qs_max = max(qs)
  print("ps_min = ", ps_min)
  print("ps_max = ", ps_max)
  print("qs_min = ", qs_min)
  print("qs_max = ", qs_max)

  analytical_times = [0.0, 1.0]
  pc_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        pc_list.append(pcs[ii])

  #print("pc = ", pc_list)
  pc_list_new = list(sorted(set(pc_list)))

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  ##Plot a
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  material_dict = get_yield_surface(uda_path)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  eqShear_vs_meanStress(ps, qs)

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('CamClay 04:\nCurved Yield Surface')
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  #Add Analytical
  #plt.legend()
  savePNG(save_path + '/Test04_verificationPlot', '1280x960')
  plt.show()


def test05_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 05 - Hydrostatic Compression Fixed Cap")
  times, sigmas = get_pTensor(uda_path, "p.stress")
  times, pcs = get_pScalar(uda_path, "p.p_c")
  ps_unscaled, qs_unscaled = get_ps_and_qs(sigmas)

  # Scale the data
  ps = [-p*1.0e-6 for p in ps_unscaled]
  qs = [-q*1.0e-6 for q in qs_unscaled]

  # Find min/max values
  ps_min = min(ps)
  ps_max = max(ps)
  qs_min = min(qs)
  qs_max = max(qs)
  print("ps_min = ", ps_min)
  print("ps_max = ", ps_max)
  print("qs_min = ", qs_min)
  print("qs_max = ", qs_max)

  analytical_times = [0.0, 1.0]
  pc_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        pc_list.append(pcs[ii])

  #print("pc = ", pc_list)
  pc_list_new = list(sorted(set(pc_list)))

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  ##Plot a
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  material_dict = get_yield_surface(uda_path)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  eqShear_vs_meanStress(ps, qs)

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('CamClay 05:\nHydrostatic Compression Fixed Cap')
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  #Add Analytical
  #plt.legend()
  savePNG(save_path + '/Test05_verificationPlot', '1280x960')
  plt.show()


def test06_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 06 - Uniaxial Strain Cap Evolution")
  times, sigmas = get_pTensor(uda_path, "p.stress")
  times, pcs = get_pScalar(uda_path, "p.p_c")
  ps_unscaled, qs_unscaled = get_ps_and_qs(sigmas)

  # Scale the data
  ps = [-p*1.0e-6 for p in ps_unscaled]
  qs = [-q*1.0e-6 for q in qs_unscaled]

  # Find min/max values
  ps_min = min(ps)
  ps_max = max(ps)
  qs_min = min(qs)
  qs_max = max(qs)
  print("ps_min = ", ps_min)
  print("ps_max = ", ps_max)
  print("qs_min = ", qs_min)
  print("qs_max = ", qs_max)

  analytical_times = [0.001, 0.2, 0.4, 0.6, 0.8, 0.999]
  pc_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        pc_list.append(pcs[ii])

  #print("pc = ", pc_list)
  pc_list_new = list(sorted(set(pc_list)))

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  ##Plot a
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  material_dict = get_yield_surface(uda_path)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  eqShear_vs_meanStress(ps, qs)

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('CamClay 06:\nUniaxial Strain Cap Evolution')
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  #Add Analytical
  #plt.legend()
  savePNG(save_path + '/Test06_verificationPlot', '1280x960')
  plt.show()


def test07_postProc(uda_path, save_path, **kwargs):

  #Extract stress history
  print("Post Processing Test: 07 - Hydrostatic Compression Cap Evolution")
  times, sigmas = get_pTensor(uda_path, "p.stress")
  times, pcs = get_pScalar(uda_path, "p.p_c")

  # Scale the data
  ps_unscaled, qs_unscaled = get_ps_and_qs(sigmas)
  ps = [-p*1.0e-6 for p in ps_unscaled]
  qs = [-q*1.0e-6 for q in qs_unscaled]

  # Find min/max values
  ps_min = min(ps)
  ps_max = max(ps)
  qs_min = min(qs)
  qs_max = max(qs)
  print("ps_min = ", ps_min)
  print("ps_max = ", ps_max)
  print("qs_min = ", qs_min)
  print("qs_max = ", qs_max)

  analytical_times = times

  # Get the material properties
  material_dict = get_yield_surface(uda_path)

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  ##Plot a
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')

  ax1 = eqShear_vs_meanStress(ps, qs)

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('CamClay 07:\nHydrostatic Compression with Fixed Cap')
  plt.ylabel(str_to_mathbf('Porosity'))
  plt.xlabel(str_to_mathbf('I_{1}:first invariant of stress tensor (MPa)'))

  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  plt.legend()
  savePNG(save_path + '/Test07_verificationPlot', '1280x960')
  plt.show()


def test08_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 08 - Loading/Unloading")
  times, sigmas = get_pTensor(uda_path, "p.stress")
  times, totalStrain = get_pTensor(uda_path, "p.strain")
  times, elasticStrain = get_pTensor(uda_path, "p.elasticStrain")
  times, pcs = get_pScalar(uda_path, "p.p_c")

  # Scale the stress data
  ps_unscaled, qs_unscaled = get_ps_and_qs(sigmas)
  ps = [-p*1.0e-6 for p in ps_unscaled]
  qs = [-q*1.0e-6 for q in qs_unscaled]
  I1s = [-3.0*p*1.0e6 for p in ps_unscaled]

  # Find min/max values
  I1s_min = min(I1s)
  I1s_max = max(I1s)
  ps_min = min(ps)
  ps_max = max(ps)
  print("I1s_min = ", I1s_min)
  print("I1s_max = ", I1s_max)
  print("ps_min = ", ps_min)
  print("ps_max = ", ps_max)

  # Compute volumetric strains
  plasticStrain = [eps - eps_e for eps in totalStrain for eps_e in elasticStrain]
  plasticStrainVol = [eps_p.trace() for eps_p in plasticStrain]
  #elasticStrainVol = [eps_e.trace() for eps_e in elasticStrain]
  totalStrainVol = [eps.trace() for eps in totalStrain]

  analytical_times = [0.0, 1.0, 2.0, 3.0, 4.0]

  ev_p_list = []
  pc_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 1.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        ev_p_list.append(plasticStrainVol[ii])
        pc_list.append(pcs[ii])

  #print("ev_p = ", ev_p_list)
  #print("pc = ", pc_list)
  #print(sorted(set(ev_p_list)))
  ev_p_list_new = list(sorted(set(ev_p_list)))

  material_dict = get_yield_surface(uda_path)

  ###PLOTTING
  int_formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  exp_formatter = ticker.FuncFormatter(exp_fmt)
  ##Plot a
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75, left=0.15)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  #print(len(times))
  #print(len(ps))
  ax1 = eqShear_vs_meanStress(times, np.array(ps))

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('CamClay 08:\nLoading/Unloading (plot a)')
  plt.ylabel(str_to_mathbf('Pressure (MPa)'))
  plt.xlabel(str_to_mathbf('Time (s)'))
  ax1.xaxis.set_major_formatter(int_formatter)
  ax1.yaxis.set_major_formatter(exp_formatter)
  ax1.tick_params(axis='both', labelsize='small')
  savePNG(save_path + '/Test08_verificationPlot_a', '1280x960')

  ##Plot b
  plt.figure(2)
  plt.clf()
  ax2 = plt.subplot(111)
  plt.subplots_adjust(right=0.75, left=0.15)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  ax1 = eqShear_vs_meanStress(times, totalStrainVol)

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('CamClay 08:\nLoading/Unloading (plot b)')
  plt.ylabel(str_to_mathbf('Total Volumetric Strain, \epsilon_{v}'))
  plt.xlabel(str_to_mathbf('Time (s)'))
  ax2.xaxis.set_major_formatter(int_formatter)
  ax2.yaxis.set_major_formatter(int_formatter)
  ax2.tick_params(axis='both', labelsize='small')
  savePNG(save_path + '/Test08_verificationPlot_b', '1280x960')

  plt.show()

#-----------------------------------------------------------------------------------
# Read the experimental stress data and compute p,q
#-----------------------------------------------------------------------------------
def readExptStressData(uda_path, file_name_expt, **kwargs):

  # Import the csv module
  import csv

  # Open the experimental data file
  I1_expt = []
  J2_expt = []
  time_expt = []
  sigma_a_expt = []
  sigma_r_expt = []
  with open(file_name_expt, 'rb') as exptfile:
    reader = csv.DictReader(exptfile, delimiter=' ', quotechar='"')
    for row in reader:
      #print(row['Axial stress (computed)'], row['mean stress (computed)'])
      time = float(row['Time'])
      sigma_a = float(row['Axial stress (computed)'])
      sigma_m = float(row['mean stress (computed)'])
      #print("sig_a = ", sigma_a, " sig_m = ", sigma_m)
      #print(type(sigma_a).__name__, type(sigma_m).__name__)

      sigma_r = 0.5 * (3.0 * sigma_m - sigma_a)
      I1 = 3.0 * sigma_m
      J2 = pow(sigma_a - sigma_r, 2) / 3.0
      #print("I1 = ", I1, " J2 = ", J2 )

      time_expt.append(time)
      sigma_a_expt.append(sigma_a)
      sigma_r_expt.append(sigma_r)
      I1_expt.append(I1)
      J2_expt.append(J2)

  pp_expt = -np.array(I1_expt) / 3.0
  qq_expt = -np.sqrt(3.0 * np.array(J2_expt))

  return time_expt, sigma_a_expt, sigma_r_expt, pp_expt, qq_expt


#-----------------------------------------------------------------------------------
# Read the simulation stress data and compute p,q (converts to MPa)
#-----------------------------------------------------------------------------------
def readSimStressData(uda_path):

  NAN_FAIL = False

  #Extract stress history
  print("Extracting stress history...")
  args = [partextract_exe, "-partvar", "p.stress", uda_path]
  print(args)
  F_stress = tempfile.TemporaryFile()
  tmp = sub_proc.Popen(args, stdout=F_stress, stderr=sub_proc.PIPE)
  dummy = tmp.wait()
  print('Done.')

  #Read file back in
  F_stress.seek(0)
  time_sim = []
  sigma_sim = []
  sigma_a_sim = []
  sigma_r_sim = []
  sigma_ar_sim = []
  for line in F_stress:

    # If the first word in the line is Error then exit
    line = line.decode()
    print("Line = ", line)
    first_word = line.partition(' ')[0]
    if first_word == "Error":
      print("**ERROR** partextract failed to read the stress history data.")
      print("          It generated the following message:")
      print(line)
      sys.exit("Stopping program.")

    line = line.strip().split()
    S11 = np.float64(line[4]) * 1.0e-6
    S12 = np.float64(line[5]) * 1.0e-6
    S13 = np.float64(line[6]) * 1.0e-6
    S21 = np.float64(line[7]) * 1.0e-6
    S22 = np.float64(line[8]) * 1.0e-6
    S23 = np.float64(line[9]) * 1.0e-6
    S31 = np.float64(line[10]) * 1.0e-6
    S32 = np.float64(line[11]) * 1.0e-6
    S33 = np.float64(line[12]) * 1.0e-6

    sigma_a = -S11
    sigma_r = -S22
    sigma_ar = -S12

    time_sim.append(float(line[0]))
    sigma_sim.append(
        np.array([[S11, S12, S13], [S21, S22, S23], [S31, S32, S33]]))
    sigma_a_sim.append(sigma_a)
    sigma_r_sim.append(sigma_r)
    sigma_ar_sim.append(sigma_ar)
    for i in range(3):
      for j in range(3):
        if np.isnan(sigma_sim[-1][i][j]):
          NAN_FAIL = True
  F_stress.close()
  if NAN_FAIL:
    print("\nERROR: 'nan's found reading in stress. Will not plot correctly")

  # Compute I1 and J2
  I1_sim = []
  J2_sim = []
  for sigma in sigma_sim:
    I1_sim.append(sigma_I1(sigma))
    J2_sim.append(sigma_J2(sigma))

  # Compute p, q
  pp_sim = np.array(I1_sim) / 3.0
  qq_sim = -np.sqrt(3.0 * np.array(J2_sim))

  return time_sim, sigma_sim, sigma_a_sim, sigma_r_sim, sigma_ar_sim, pp_sim, qq_sim


#---------------------------------------------------------------------------------
# Interpolate the data at a set of given time points
#---------------------------------------------------------------------------------
def getDataTimeSnapshots(time_snapshots, time_expt, data_expt):

  # Create a clean list containing the data as functions of time
  time_list = []
  val_list = []
  for ii in range(1, len(time_expt)):
    tt_0 = time_expt[ii - 1]
    tt_1 = time_expt[ii]
    for jj, ta in enumerate(time_snapshots):
      ss = (ta - tt_0) / (tt_1 - tt_0)
      if (ss >= 0.0 and ss < 1.0):
        #print("ii = " , ii, " tt_0 = " , tt_0, " tt_1 = ", tt_1, " jj = " , jj, " ta = " , ta )
        time_list.append(ta)
        val_list.append((1 - ss) * data_expt[ii - 1] + ss * data_expt[ii])

  return time_list, val_list


#-----------------------------------------------------------------------------------
# Plot the expt data as a function of time (aasume MPa)
#-----------------------------------------------------------------------------------
def plotExptDataSigmaTime(fig, time_snapshots, time_expt, sigma_a_expt,
                          sigma_r_expt):

  # Get snapshots from expt data
  time_snap, sigma_a_snap = getDataTimeSnapshots(time_snapshots, time_expt,
                                                 sigma_a_expt)
  time_snap, sigma_r_snap = getDataTimeSnapshots(time_snapshots, time_expt,
                                                 sigma_r_expt)

  # Activate the figure
  plt.figure(fig.number)

  # Plot sigma_a vs. time
  plt.plot(time_expt, sigma_a_expt, '-r', label='$\sigma_a$ (expt)')
  plt.plot(time_expt, sigma_r_expt, '-b', label='$\sigma_r$ (expt)')

  # Plot filled circles at time snapshots
  for ii in range(0, len(time_snap)):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(time_snap))
    #xvals = np.array((time_snap[ii], time_snap[ii+1]))
    #yvals_a = np.array((sigma_a_snap[ii], sigma_a_snap[ii+1]))
    #yvals_r = np.array((sigma_r_snap[ii], sigma_r_snap[ii+1]))
    #plt.plot(xvals, yvals_a, '-', color=plt_color)
    #plt.plot(xvals, yvals_r, '-', color=plt_color)
    plt.plot(time_snap[ii], sigma_a_snap[ii], 'v', color=plt_color)
    plt.plot(time_snap[ii], sigma_r_snap[ii], 'v', color=plt_color)

  return time_snap, sigma_a_snap, sigma_r_snap


#-----------------------------------------------------------------------------------
# Plot the sim data as a function of time (assume MPa)
#-----------------------------------------------------------------------------------
def plotSimDataSigmaTime(fig, time_snapshots, time_sim, sigma_a_sim,
                         sigma_r_sim, sigma_ar_sim):

  # Get snapshots from expt data
  time_snap, sigma_a_snap = getDataTimeSnapshots(time_snapshots, time_sim,
                                                 sigma_a_sim)
  time_snap, sigma_r_snap = getDataTimeSnapshots(time_snapshots, time_sim,
                                                 sigma_r_sim)
  time_snap, sigma_ar_snap = getDataTimeSnapshots(time_snapshots, time_sim,
                                                  sigma_ar_sim)

  # Activate the figure
  plt.figure(fig.number)

  # Plot sigma_a vs. time
  plt.plot(time_sim, sigma_a_sim, '--r', label='$\sigma_a$ (sim)')
  plt.plot(time_sim, sigma_r_sim, '--b', label='$\sigma_r$ (sim)')
  plt.plot(time_sim, sigma_ar_sim, '--g', label='$\sigma_{ar}$ (sim)')

  # Plot filled circles at time snapshots
  for ii in range(0, len(time_snap)):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(time_snap))
    plt.plot(time_snap[ii], sigma_a_snap[ii], 'o', color=plt_color)
    plt.plot(time_snap[ii], sigma_r_snap[ii], 'o', color=plt_color)
    #plt.plot(time_snap[ii], sigma_ar_snap[ii], 'o', color=plt_color)

  return time_snap, sigma_a_snap, sigma_r_snap, sigma_ar_snap


#-----------------------------------------------------------------------------------
# Plot the expt data (pq) as a function of time (aasume MPa)
#-----------------------------------------------------------------------------------
def plotExptDataPQTime(fig, time_snapshots, time_expt, p_expt, q_expt):

  # Get snapshots from expt data
  time_snap, p_snap = getDataTimeSnapshots(time_snapshots, time_expt, p_expt)
  time_snap, q_snap = getDataTimeSnapshots(time_snapshots, time_expt, q_expt)

  # Activate the figure
  plt.figure(fig.number)

  # Plot sigma_a vs. time
  plt.plot(time_expt, p_expt, '-r', label='$p$ (expt)')
  plt.plot(time_expt, q_expt, '-b', label='$q$ (expt)')

  # Plot filled circles at time snapshots
  for ii in range(0, len(time_snap)):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(time_snap))
    plt.plot(time_snap[ii], p_snap[ii], 'v', color=plt_color)
    plt.plot(time_snap[ii], q_snap[ii], 'v', color=plt_color)

  return time_snap, p_snap, q_snap


#-----------------------------------------------------------------------------------
# Plot the sim data (pq) as a function of time (aasume MPa)
#-----------------------------------------------------------------------------------
def plotSimDataPQTime(fig, time_snapshots, time_sim, p_sim, q_sim):

  # Get snapshots from sim data
  time_snap, p_snap = getDataTimeSnapshots(time_snapshots, time_sim, p_sim)
  time_snap, q_snap = getDataTimeSnapshots(time_snapshots, time_sim, q_sim)

  # Activate the figure
  plt.figure(fig.number)

  # Plot sigma_a vs. time
  plt.plot(time_sim, p_sim, '--r', label='$p$ (sim)')
  plt.plot(time_sim, q_sim, '--b', label='$q$ (sim)')

  # Plot filled circles at time snapshots
  for ii in range(0, len(time_snap)):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(time_snap))
    plt.plot(time_snap[ii], p_snap[ii], 'o', color=plt_color)
    plt.plot(time_snap[ii], q_snap[ii], 'o', color=plt_color)

  return time_snap, p_snap, q_snap


#---------------------------------------------------------------------------------
# Read the internal state variables from the UDA file:
#    1) volumetric plastic strain
#    2) capX
#    3) kappa
#    3) zeta
#---------------------------------------------------------------------------------
def getAllInternalVariables(uda_path):

  # Get the internal variables
  time_list, pc_list = get_pScalar(uda_path, "p.p_c")

  return pc_list, time_list


#---------------------------------------------------------------------------------
# Read the internal state variables from the UDA file:
#    1) elastic strain
#    2) total strain
#    2) pc
#---------------------------------------------------------------------------------
def getInternalVariables(uda_path, analytical_times):

  # Get the internal variables
  pc, times = getAllInternalVariables(uda_path)

  # Create a clean list containing the data as functions of time
  pc_list = []
  time_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  #print("Analytical times = ", an_times_add)
  for ii in range(1, len(times)):
    tt_0 = times[ii - 1]
    tt_1 = times[ii]
    for jj, ta in enumerate(an_times_add):
      ss = (ta - tt_0) / (tt_1 - tt_0)
      #print("ii = " , ii, " tt_0 = " , tt_0, " tt_1 = ", tt_1, " jj = " , jj, " ta = " , ta, " ss = ", ss )
      if (ss >= 0.0 and ss < 1.0):
        #print("ii = " , ii, " tt_0 = " , tt_0, " tt_1 = ", tt_1, " jj = " , jj, " ta = " , ta )
        time_list.append(ta)
        pc_list.append((1 - ss) * pc[ii - 1] + ss * pc[ii])
    #    break
    #else:
    #  continue

  #print("time = ", time_list)
  #print("pc = ", pc_list)

  return pc_list, time_list


#-----------------------------------------------------------------------------
# Plot the yield surface based on the compute internal variables
#-----------------------------------------------------------------------------
def plotPQYieldSurfaceSim(uda_path, time_points, **kwargs):

  # Read the keyword arguments
  p_min_expt = 0
  q_min_expt = 0
  if 'p_min_expt' in kwargs:
    p_min_expt = kwargs['p_min_expt']
  if 'q_min_expt' in kwargs:
    q_min_expt = kwargs['q_min_expt']

  # get the internal variables
  pc_list, time_list = getInternalVariables(uda_path, time_points)

  # Get material parameters
  material_dict = get_yield_surface(uda_path)
  p0 = material_dict['p0']
  alpha = material_dict['alpha']
  kappatilde = material_dict['kappatilde']
  epse_v0 = material_dict['epse_v0']
  mu0 = material_dict['mu0']
  M = material_dict['M']
  pc0 = material_dict['pc0']
  lambdatilde = material_dict['lambdatilde']

  # Find pmin and qmin
  pmin = 1.0e8
  qmax = -1.0e8

  # Loop through time snapshots
  for ii, pc in enumerate(pc_list):

    # Get the internal variables
    pc = pc_list[ii]

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(pc_list))

    # Create an array of I1 values
    num_points = 100
    ps = np.linspace(p0, pc, num_points)
    qs = []

    #q versus p
    for p in ps:

      # Compute q
      q = M * np.sqrt(np.abs(p * (p - pc)))
      qs.append(q)

      #print("q = ", q, "p = ", p, "pc = ", pc)

    xs = [-p * 1.0e-6 for p in ps]
    ys = [q * 1.0e-6 for q in qs]
    xs = np.array(xs)
    ys = np.array(ys)

    if (min(xs) < pmin):
      pmin = min(xs)

    if (max(ys) > qmax):
      qmax = max(ys)

    #print(xs.shape, ys.shape)
    pc_str = "%.2g" % (pc)
    #print("pc " , pc_str)
    #time_str = "%.2g" % time_list[ii]
    #print("time " , time_str)

    label_str = ' $p_c$ = ' + pc_str
    line1 = plt.plot(xs, ys, '--b', linewidth=1, label=label_str)
    line2 = plt.plot(xs, -ys, '--b', linewidth=1)
    plt.setp(line1, color=plt_color)
    plt.setp(line2, color=plt_color)
    plt.legend(loc=2, prop={'size': 8})

  print("pmin = ", pmin, " qmax = ", qmax)
  axes = plt.gca()
  #axes.set_xlim([1.3*pmin, 1])
  #axes.set_ylim([-1.3*qmax, 1.3*qmax])
  #axes.set_xlim([-1, -1.3 * pmin])
  #axes.set_ylim([-1.3 * qmax, 1.3 * qmax])
  return pmin, qmax

