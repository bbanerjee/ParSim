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


def get_pStress(uda_path):
  NAN_FAIL = False
  #Extract stress history
  print("Extracting stress history...")
  args = [partextract_exe, "-partvar", "p.stress", uda_path]
  print(args)
  F_stress = tempfile.TemporaryFile()
  #F_stress = open("./tempStressFileOut.txt","w+")
  #open(os.path.split(uda_path)[0]+'/stressHistory.dat',"w+")
  tmp = sub_proc.Popen(args, stdout=F_stress, stderr=sub_proc.PIPE)
  dummy = tmp.wait()
  print('Done.')
  #Read file back in
  F_stress.seek(0)
  times = []
  sigmas = []
  for line in F_stress:

    # If the first word in the line is Error then exit
    #print("Line = ", line)
    first_word = line.partition(' ')[0]
    if first_word == "Error":
      print("**ERROR** partextract failed to read the stress history data.")
      print("          It generated the following message:")
      print(line)
      sys.exit("Stopping program.")

    line = line.strip().split()
    times.append(float(line[0]))
    S11 = np.float64(line[4])
    S12 = np.float64(line[5])
    S13 = np.float64(line[6])
    S21 = np.float64(line[7])
    S22 = np.float64(line[8])
    S23 = np.float64(line[9])
    S31 = np.float64(line[10])
    S32 = np.float64(line[11])
    S33 = np.float64(line[12])
    sigmas.append(np.array([[S11, S12, S13], [S21, S22, S23], [S31, S32, S33]]))
    for i in range(3):
      for j in range(3):
        if np.isnan(sigmas[-1][i][j]):
          NAN_FAIL = True
  F_stress.close()
  if NAN_FAIL:
    print("\nERROR: 'nan's found reading in stress. Will not plot correctly")
  return times, sigmas


def get_pDeformationMeasure(uda_path):
  NAN_FAIL = False
  #Extract stress history
  print("Extracting deformation history...")
  args = [partextract_exe, "-partvar", "p.deformationGradient", uda_path]
  F_defMes = tempfile.TemporaryFile()
  #open(os.path.split(uda_path)[0]+'/stressHistory.dat',"w+")
  tmp = sub_proc.Popen(args, stdout=F_defMes, stderr=sub_proc.PIPE)
  dummy = tmp.wait()
  print('Done.')
  #Read file back in
  F_defMes.seek(0)
  times = []
  Fs = []
  for line in F_defMes:
    line = line.strip().split()
    times.append(float(line[0]))
    F11 = np.float64(line[4])
    F12 = np.float64(line[5])
    F13 = np.float64(line[6])
    F21 = np.float64(line[7])
    F22 = np.float64(line[8])
    F23 = np.float64(line[9])
    F31 = np.float64(line[10])
    F32 = np.float64(line[11])
    F33 = np.float64(line[12])
    Fs.append(np.array([[F11, F12, F13], [F21, F22, F23], [F31, F32, F33]]))
    for i in range(3):
      for j in range(3):
        if np.isnan(Fs[-1][i][j]):
          NAN_FAIL = True
  F_defMes.close()
  if NAN_FAIL:
    print("\nERROR: 'nan's found reading in stress. Will not plot correctly")
  return times, Fs


def get_epsilons(uda_path):
  #Assumes no shear strains
  times, Fs = get_pDeformationMeasure(uda_path)
  epsils = []
  for F in Fs:
    epsils.append(
        np.array([[np.log(F[0][0]), 0, 0], [0, np.log(F[1][1]), 0],
                  [0, 0, np.log(F[2][2])]]))
  return times, epsils


def get_pKappa(uda_path):
  #Extract stress history
  print("Extracting kappa history...")
  args = [partextract_exe, "-partvar", "p.kappa", uda_path]
  F_kappa = tempfile.TemporaryFile()
  #open(os.path.split(uda_path)[0]+'/kappaHistory.dat',"w+")
  tmp = sub_proc.Popen(args, stdout=F_kappa, stderr=sub_proc.PIPE)
  dummy = tmp.wait()
  print('Done.')
  #Read file back in
  F_kappa.seek(0)
  times = []
  kappas = []
  for line in F_kappa:
    line = line.strip().split()
    times.append(float(line[0]))
    kappas.append(float(line[4]))
  F_kappa.close()
  return times, kappas


def get_pPlasticStrainVol(uda_path):
  FAIL_NAN = False
  #Extract stress history
  print("Extracting plasticStrainVol history...")
  args = [partextract_exe, "-partvar", "p.evp", uda_path]
  F_plasticStrainVol = tempfile.TemporaryFile()
  #open(os.path.split(uda_path)[0]+'/plasticStrainVolHistory.dat',"w+")
  tmp = sub_proc.Popen(args, stdout=F_plasticStrainVol, stderr=sub_proc.PIPE)
  dummy = tmp.wait()
  print('Done.')
  #Read file back in
  F_plasticStrainVol.seek(0)
  times = []
  plasticStrainVol = []
  for line in F_plasticStrainVol:
    line = line.strip().split()
    times.append(float(line[0]))
    plasticStrainVol.append(np.float64(line[4]))
    if np.isnan(plasticStrainVol[-1]):
      FAIL_NAN = True
  F_plasticStrainVol.close()
  if FAIL_NAN:
    print(
        "\ERROR: 'nan' encountered while retrieving p.evp, will not plot correctly."
    )
  return times, plasticStrainVol


def get_pElasticStrainVol(uda_path):
  FAIL_NAN = False
  #Extract elastic strain history
  print("Extracting elasticStrainVol history...")
  args = [partextract_exe, "-partvar", "p.eve", uda_path]
  F_elasticStrainVol = tempfile.TemporaryFile()
  #open(os.path.split(uda_path)[0]+'/elasticStrainVolHistory.dat',"w+")
  tmp = sub_proc.Popen(args, stdout=F_elasticStrainVol, stderr=sub_proc.PIPE)
  dummy = tmp.wait()
  print('Done.')
  #Read file back in
  F_elasticStrainVol.seek(0)
  times = []
  elasticStrainVol = []
  for line in F_elasticStrainVol:
    line = line.strip().split()
    times.append(float(line[0]))
    elasticStrainVol.append(np.float64(line[4]))
    if np.isnan(elasticStrainVol[-1]):
      FAIL_NAN = True
  F_elasticStrainVol.close()
  if FAIL_NAN:
    print(
        "\ERROR: 'nan' encountered while retrieving p.eve, will not plot correctly."
    )
  return times, elasticStrainVol


def get_totalStrainVol(uda_path):
  times, plasticStrainVol = get_pPlasticStrainVol(uda_path)
  times, elasticStrainVol = get_pElasticStrainVol(uda_path)

  print('num plastic : ', len(plasticStrainVol))
  print('num elastic : ', len(elasticStrainVol))

  totalStrainVol = np.array(plasticStrainVol) + np.array(elasticStrainVol)
  return times, totalStrainVol


def get_capX(uda_path):
  #Extract capX history
  print("Extracting capX history...")
  args = [partextract_exe, "-partvar", "p.CapX", uda_path]
  F_capX = tempfile.TemporaryFile()
  #open(os.path.split(uda_path)[0]+'/kappaHistory.dat',"w+")
  tmp = sub_proc.Popen(args, stdout=F_capX, stderr=sub_proc.PIPE)
  dummy = tmp.wait()
  print('Done.')
  #Read file back in
  F_capX.seek(0)
  times = []
  capX = []
  for line in F_capX:
    line = line.strip().split()
    times.append(float(line[0]))
    capX.append(float(line[4]))
  F_capX.close()
  return times, capX


def get_zeta(uda_path):
  #Extract zeta history
  print("Extracting zeta history...")
  args = [partextract_exe, "-partvar", "p.Zeta", uda_path]
  F_zeta = tempfile.TemporaryFile()
  #open(os.path.split(uda_path)[0]+'/kappaHistory.dat',"w+")
  tmp = sub_proc.Popen(args, stdout=F_zeta, stderr=sub_proc.PIPE)
  dummy = tmp.wait()
  print('Done.')
  #Read file back in
  F_zeta.seek(0)
  times = []
  zeta = []
  for line in F_zeta:
    line = line.strip().split()
    times.append(float(line[0]))
    zeta.append(float(line[4]))
  F_zeta.close()
  return times, zeta


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
                          LINE_LABEL='Uintah',
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
  #Reads in FSLOPE, FSLOPE_p, PEAKI1, CR, and P0
  #WILL ONLY WORK FOR SINGLE ELEMENT TESTS OR DECKS
  #HAVING ONLY ONE ARENISCA SPECIFICATION
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
  print(sorted_list_tuples)
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


def get_kappa(PEAKI1, P0, CR):
  PEAKI1, P0, CR
  kappa = PEAKI1 - CR * (PEAKI1 - P0)
  return kappa


def get_rs(nPoints, FSLOPE, PEAKI1, P0, CR):

  kappa = get_kappa(PEAKI1, P0, CR)
  I1s = np.linspace(PEAKI1, P0, nPoints)
  rs = []
  for I1 in I1s:
    inner_root = (1.0 - (pow(kappa - I1, 2.0) / pow(kappa - P0, 2.0)))
    r = FSLOPE * (I1 - PEAKI1) * np.sqrt(2.0 * inner_root)
    rs.append(r)
  return I1s, rs


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


def plot_crush_curve(uda_path, I1lims=[-10000, 0]):
  nPoints = 500
  material_dict = get_yield_surface(uda_path)
  P0 = material_dict['P0']
  P1 = material_dict['P1']
  P3 = material_dict['P3']

  # Analytical solution for porosity vs. X for piece-wise crush curve
  # compression
  I1sC = np.linspace(I1lims[0] * 1.0e6, P0, nPoints)
  porosityC = 1 - np.exp(-P3 * np.exp(P1 * (I1sC - P0)))
  # tension
  I1sT = np.linspace(P0, I1lims[1] * 1.0e6, nPoints)
  porosityT = 1 - np.exp(-(I1sT / P0)**(P0 * P1 * P3) - P3 + 1)

  # Scale down again
  I1_c = []
  for ii in I1sC:
    I1_c.append(ii * 1.0e-6)
  I1_t = []
  for ii in I1sT:
    I1_t.append(ii * 1.0e-6)

  plt.plot(I1_c,
           porosityC,
           '--g',
           linewidth=lineWidth + 1,
           label='Analytical crush curve - Compression')
  plt.hold(True)
  plt.plot(I1_t,
           porosityT,
           '--b',
           linewidth=lineWidth + 1,
           label='Analytical crush curve - Tension')
  #plt.ylabel(str_to_mathbf('I_{1} (MPa)'))
  #plt.xlabel(str_to_mathbf('Porosity'))


def plot_yield_surface_OLD(uda_path):
  nPoints = 500
  material_dict = get_yield_surface(uda_path)
  FSLOPE = material_dict['FSLOPE']
  #FSLOPE_p = material_dict['FSLOPE']
  PEAKI1 = material_dict['PEAKI1']
  CR = material_dict['CR']
  P0 = material_dict['P0']
  I1s, rs = get_rs(nPoints, FSLOPE, PEAKI1, P0, CR)
  zbars = I1_to_zbar(I1s)
  #WTF?
  for i in range(len(rs)):
    rs[i] = -rs[i]
  #print(zbars)
  #print(rs)
  plt.plot(np.array(I1s) / 3.0,
           rs,
           '--k',
           linewidth=lineWidth + 1,
           label='Initial Yield Surface')
  plt.plot(np.array(I1s) / 3.0, -np.array(rs), '--k', linewidth=lineWidth + 1)


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


def J2_at_Yield(uda_path):
  material_dict = get_yield_surface(uda_path)
  PEAKI1 = material_dict['PEAKI1']
  FSLOPE = material_dict['FSLOPE']
  #FSLOPE_p = material_dict['FSLOPE']
  STREN = material_dict['STREN']
  YSLOPE = material_dict['YSLOPE']
  BETA = material_dict['BETA']
  B0 = material_dict['B0']
  B1 = material_dict['B1']
  B2 = material_dict['B2']
  B3 = material_dict['B3']
  B4 = material_dict['B4']
  G0 = material_dict['G0']
  P0 = material_dict['P0']
  P1 = material_dict['P1']
  P3 = material_dict['P3']
  CR = material_dict['CR']
  fluid_B0 = material_dict['fluid_B0']
  Pf0 = material_dict['P_f0']
  T1 = material_dict['T1']
  T2 = material_dict['T2']
  subcyc_char_num = material_dict['subcycling char num']
  #hardening_const = material_dict['hardening_constant']

  kappa_initial = get_kappa(PEAKI1, P0, CR)
  I1 = 0
  I1_plus3Pf0 = I1 + 3.0 * Pf0
  if I1_plus3Pf0 >= kappa_initial and I1_plus3Pf0 <= PEAKI1:
    J2 = (FSLOPE**2) * ((I1 - PEAKI1 + 3.0 * Pf0)**2)
  elif I1_plus3Pf0 >= P0 and I1_plus3Pf0 < kappa_initial:
    J2 = ((FSLOPE**2) * ((I1 - PEAKI1 + 3.0 * Pf0)**2)) * (1.0 - (
        (I1 + CR * FSLOPE * I1 - P0 - CR * FSLOPE * PEAKI1 + 3.0 * Pf0 +
         3.0 * CR * FSLOPE * Pf0)**2 / ((CR**2) * (FSLOPE**2) *
                                        (P0 - PEAKI1)**2)))
  else:
    J2 = 0.0
  return J2


def plot_yield_surface(uda_path, PLOT_TYPE='J2_vs_I1'):
  num_points = 500
  material_dict = get_yield_surface(uda_path)
  PEAKI1 = material_dict['PEAKI1']
  FSLOPE = material_dict['FSLOPE']
  #FSLOPE_p = material_dict['FSLOPE']
  STREN = material_dict['STREN']
  YSLOPE = material_dict['YSLOPE']
  BETA = material_dict['BETA']
  B0 = material_dict['B0']
  B1 = material_dict['B1']
  B2 = material_dict['B2']
  B3 = material_dict['B3']
  B4 = material_dict['B4']
  G0 = material_dict['G0']
  P0 = material_dict['P0']
  P1 = material_dict['P1']
  P3 = material_dict['P3']
  CR = material_dict['CR']
  fluid_B0 = material_dict['fluid_B0']
  Pf0 = material_dict['P_f0']
  T1 = material_dict['T1']
  T2 = material_dict['T2']
  subcyc_char_num = material_dict['subcycling char num']
  #hardening_const = material_dict['hardening_constant']
  kappa_initial = get_kappa(PEAKI1, P0, CR)
  I1s = np.linspace(P0 - 3.0 * Pf0, PEAKI1 - 3.0 * Pf0, num_points)
  #print('Region 1:: ','I1 >= kappa initial-3.0*Pf0 : ',kappa_initial-3.0*Pf0,' ','I1 <= PEAKI1-3*Pf0 : ',PEAKI1-3.0*Pf0)
  #print('Region 2:: ','I1 >= P0-3*Pf0 : ',P0-3.0*Pf0,' ','I1 < kappa_initial-3*Pf0 : ',kappa_initial-3.0*Pf0)
  #print('Region 3:: Not Region 1 or 2')

  #J2 versus I1
  J2s = []
  PLOT = True
  for I1 in I1s:
    I1_plus3Pf0 = I1 + 3.0 * Pf0
    if I1_plus3Pf0 >= kappa_initial and I1_plus3Pf0 <= PEAKI1:
      J2 = (FSLOPE**2) * ((I1 - PEAKI1 + 3.0 * Pf0)**2)
    elif I1_plus3Pf0 >= P0 and I1_plus3Pf0 < kappa_initial:
      Ff = FSLOPE * (PEAKI1 - I1_plus3Pf0)
      fc = np.sqrt(1.0 - ((kappa_initial - I1_plus3Pf0) /
                          (kappa_initial - P0))**2)
      J2 = (Ff * fc)**2
    else:
      J2 = 0.0
    J2s.append(J2)

  if PLOT_TYPE == 'J2_vs_I1':
    xs = np.array(I1s) * 1.0e-6
    ys = np.array(J2s) * 1.0e-12
  elif PLOT_TYPE == 'sqrtJ2_vs_I1':
    xs = np.array(I1s) * 1.0e-6
    ys = np.sqrt(np.array(J2s)) * 1.0e-6
  elif PLOT_TYPE == 'r_vs_z':
    xs = np.array(I1s) / np.sqrt(3.0) * 1.0e-6
    ys = np.sqrt(2.0 * np.array(J2s)) * 1.0e-6
  elif PLOT_TYPE == 'q_vs_I1':
    xs = np.array(I1s) * 1.0e-6
    ys = np.sqrt(3.0 * np.array(J2s)) * 1.0e-6
  elif PLOT_TYPE == 'q_vs_p':
    xs = np.array(I1s) / 3.0 * 1.0e-6
    ys = np.sqrt(3.0 * np.array(J2s)) * 1.0e-6
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
  plot_yield_surface_2(uda_path, 'J2_vs_I1')
  plt.show()
  plot_yield_surface_2(uda_path, 'sqrtJ2_vs_I1')
  plt.show()
  plot_yield_surface_2(uda_path, 'r_vs_z')
  plt.show()
  plot_yield_surface_2(uda_path, 'q_vs_I1')
  plt.show()
  plot_yield_surface_2(uda_path, 'q_vs_p')
  plt.show()


def plot_yield_surface_updated(uda_path,
                               ev_p,
                               zeta=0.0,
                               PLOT_TYPE='J2_vs_I1',
                               COLOR='#0000ff'):
  num_points = 500
  material_dict = get_yield_surface(uda_path)
  PEAKI1 = material_dict['PEAKI1']
  FSLOPE = material_dict['FSLOPE']
  #FSLOPE_p = material_dict['FSLOPE']
  STREN = material_dict['STREN']
  YSLOPE = material_dict['YSLOPE']
  BETA = material_dict['BETA']
  B0 = material_dict['B0']
  B1 = material_dict['B1']
  B2 = material_dict['B2']
  B3 = material_dict['B3']
  B4 = material_dict['B4']
  G0 = material_dict['G0']
  P0 = material_dict['P0']
  P1 = material_dict['P1']
  P3 = material_dict['P3']
  CR = material_dict['CR']
  fluid_B0 = material_dict['fluid_B0']
  Pf0 = material_dict['P_f0']
  T1 = material_dict['T1']
  T2 = material_dict['T2']
  subcyc_char_num = material_dict['subcycling char num']
  #hardening_const = material_dict['hardening_constant']

  # Set up constants
  a1 = STREN
  a2 = (FSLOPE - YSLOPE) / (STREN - YSLOPE * PEAKI1)
  a3 = (STREN - YSLOPE * PEAKI1) * math.exp(-a2 * PEAKI1)
  a4 = YSLOPE
  phi_i = 1.0 - math.exp(-P3)
  Km = B0 + B1
  Kf = fluid_B0
  C1 = Kf * (1.0 - phi_i) + Km * (phi_i)
  print('Kf = ', Kf, ' Km = ', Km)
  ev0 = 0.0
  if (Kf != 0):
    ev0 = C1 * Pf0 / (Kf * Km)

  # compute capX
  cap_X = computeCapX(ev_p, P3, P0, P1, Kf, Km, ev0, C1, phi_i, B0, B1, B2, B3,
                      B4)
  print("capX = ", cap_X)

  # Compute kappa
  kappa = PEAKI1 - CR * (PEAKI1 - cap_X)
  kappa_initial = PEAKI1 - CR * (PEAKI1 - P0)
  print("kappa = ", kappa, " kappa_initial = ", kappa_initial)

  #I1s = np.linspace(P0-3.0*Pf0,PEAKI1-3.0*Pf0,num_points)
  I1s = np.linspace(0.999 * cap_X - 3.0 * Pf0, PEAKI1 - 3.0 * Pf0, num_points)
  #print('Region 1:: ','I1 >= kappa initial-3.0*Pf0 : ',kappa_initial-3.0*Pf0,' ','I1 <= PEAKI1-3*Pf0 : ',PEAKI1-3.0*Pf0)
  #print('Region 2:: ','I1 >= P0-3*Pf0 : ',P0-3.0*Pf0,' ','I1 < kappa_initial-3*Pf0 : ',kappa_initial-3.0*Pf0)
  #print('Region 3:: Not Region 1 or 2')

  #J2 versus I1
  J2s = []
  PLOT = True
  for I1 in I1s:
    I1f = I1 + 3.0 * Pf0

    # Compute F_f
    I1_minus_zeta = I1f - zeta
    Ff = a1 - a3 * math.exp(a2 * I1_minus_zeta) - a4 * (I1_minus_zeta)
    Ff_sq = Ff**2

    # Compute Fc
    Fc_sq = 1.0
    if (I1_minus_zeta < kappa) and (cap_X < I1_minus_zeta):
      ratio = (kappa - I1_minus_zeta) / (kappa - cap_X)
      Fc_sq = 1.0 - ratio**2

    # Compute J2
    J2 = Ff_sq * Fc_sq
    J2s.append(J2)

  print(len(I1s), len(J2s))
  if PLOT_TYPE == 'J2_vs_I1':
    xs = np.array(I1s) * 1.0e-6
    ys = np.array(J2s) * 1.0e-12
  elif PLOT_TYPE == 'sqrtJ2_vs_I1':
    xs = np.array(I1s) * 1.0e-6
    ys = np.sqrt(np.array(J2s)) * 1.0e-6
  elif PLOT_TYPE == 'r_vs_z':
    xs = np.array(I1s) / np.sqrt(3.0) * 1.0e-6
    ys = np.sqrt(2.0 * np.array(J2s)) * 1.0e-6
  elif PLOT_TYPE == 'q_vs_I1':
    xs = np.array(I1s) * 1.0e-6
    ys = np.sqrt(3.0 * np.array(J2s)) * 1.0e-6
  elif PLOT_TYPE == 'q_vs_p':
    xs = np.array(I1s) / 3.0 * 1.0e-6
    ys = np.sqrt(3.0 * np.array(J2s)) * 1.0e-6
  else:
    PLOT = False
    print(
        '\nError: invalid plot type specified for updated yield surface plot.\n\tPLOT_TYPE:',
        PLOT_TYPE)
  if PLOT:
    print(xs.shape, ys.shape)
    ev_p_str = str(ev_p)
    label_str = 'Yield surf. evp = ' + ev_p_str
    line1 = plt.plot(xs, ys, '--b', linewidth=lineWidth + 1, label=label_str)
    line2 = plt.plot(xs, -ys, '--b', linewidth=lineWidth + 1)
    plt.setp(line1, color=COLOR)
    plt.setp(line2, color=COLOR)


def test_yield_surface(uda_path):
  plot_yield_surface_2(uda_path, 'J2_vs_I1')
  plt.show()
  plot_yield_surface_2(uda_path, 'sqrtJ2_vs_I1')
  plt.show()
  plot_yield_surface_2(uda_path, 'r_vs_z')
  plt.show()
  plot_yield_surface_2(uda_path, 'q_vs_I1')
  plt.show()
  plot_yield_surface_2(uda_path, 'q_vs_p')
  plt.show()


###
# Compute capX
###
def computeCapX(ev_p, P3, P0, P1, Kf, Km, ev_0, C1, phi_i, B0, B1, B2, B3, B4):

  if (ev_p <= -P3):
    return P0 * 1.0e12
  if (ev_p <= 0.0):
    capX = (P0 * P1 + math.log(1.0 + ev_p / P3)) / P1
  else:
    capX = P0 * math.pow(1.0 + ev_p, 1.0 / (P0 * P1 * P3))

  if ((Kf != 0.0) and (ev_p <= ev_0)):
    K_dry = B0 + B1 * math.exp(2.0 * B2 / capX)
    #
    if (ev_p < 0.0):
      K_dry = K_dry - B3 * math.exp(B4 / ev_p)
    #
    K_sat, shear = computeElastic(0.5 * capX, ev_p, Kf, Km, ev_0, C1, phi_i, B0,
                                  B1, B2, B3, B4, 0)
    capX = capX * K_sat / K_dry

  return capX


###
# Compute elastic
###
def computeElastic(I1, ev_p, Kf, Km, ev0, C1, phi_i, B0, B1, B2, B3, B4, G0):

  bulk = B0
  shear = G0

  if (ev_p <= 0.0):

    if (ev_p < 0.0):
      bulk = bulk - B3 * math.exp(B4 / ev_p)

    if (I1 < 0.0):
      bulk = bulk + B1 * math.exp(B2 / I1)

  if (ev_p <= ev0 and Kf != 0.0):
    Kd = B0
    if (ev_p < 0.0):
      Kd = B0 - B3 * math.exp(B4 / ev_p)
      C2 = math.exp(ev_p * Km / C1) * phi_i
      phi = C2 / (-math.exp(ev_p * Kf / C1) * (phi_i - 1.0) + C2)
      tmp = 1.0 - Kd / Km
      bulk = Kd + tmp**2 / ((tmp - phi) / Km + (1.0 / Kf - 1.0 / Km) * phi)

  return bulk, shear


### ----------
#  Test Methods Below
### ----------
def test00_postProc(uda_path, save_path, expt_stress_file, **kwargs):
  print("PostProc Test: 00 - Uniaxial Compression SHPB")

  # Just get the name of the file without extension
  file_name = os.path.splitext(expt_stress_file)[0]

  # Read the experimental data
  time_expt, sigma_a_expt, sigma_r_expt, pp_expt, qq_expt = readExptStressData(
      uda_path, expt_stress_file)

  # Read the simulation data
  times, sigmas, sigma_a_sim, sigma_r_sim, sigma_ar_sim, pp_sim, qq_sim = readSimStressData(
      uda_path)

  # Get the model parameters
  material_dict = get_yield_surface(uda_path)
  param_text = material_dict['material string']

  # Set up time points
  analytical_times = np.linspace(0.0, 0.001, 10)

  # Get snapshots of pq data (both expt and sim)
  t_expt_snap, p_expt_snap = getDataTimeSnapshots(analytical_times, time_expt,
                                                  pp_expt)
  t_expt_snap, q_expt_snap = getDataTimeSnapshots(analytical_times, time_expt,
                                                  qq_expt)
  t_sim_snap, p_sim_snap = getDataTimeSnapshots(analytical_times, times, pp_sim)
  t_sim_snap, q_sim_snap = getDataTimeSnapshots(analytical_times, times, qq_sim)

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')

  #----------------------------------------------------------------
  # Plot the yield surface for test1
  #----------------------------------------------------------------
  # Set up figure
  fig2 = plt.figure(2)
  plt.clf()
  plt.subplots_adjust(right=0.75)
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')

  # Plot p vs. q simulation results
  eqShear_vs_meanStress(-pp_sim, -qq_sim)

  # Plot filled circles at time snapshots
  for ii in range(0, len(t_sim_snap)):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(t_sim_snap))
    #xvals = np.array((p_sim_snap[ii], p_sim_snap[ii+1]))
    #yvals = np.array((q_sim_snap[ii], q_sim_snap[ii+1]))
    #plt.plot(xvals, yvals, '-', color=plt_color)
    plt.plot(-p_sim_snap[ii], -q_sim_snap[ii], 'o', color=plt_color)

  # Plot the experimental data
  line1 = plt.plot(-pp_expt, -qq_expt, '--b', linewidth=2, label='Expt.')

  # Plot filled circles at time snapshots
  for ii in range(0, len(t_expt_snap)):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(t_expt_snap))
    #xvals = np.array((p_expt_snap[ii], p_expt_snap[ii+1]))
    #yvals = np.array((q_expt_snap[ii], q_expt_snap[ii+1]))
    #plt.plot(xvals, yvals, '-', color=plt_color)
    plt.plot(-p_expt_snap[ii], -q_expt_snap[ii], 'v', color=plt_color)

  # Plot yield surfaces
  pp_expt_min = min(pp_expt)
  qq_expt_min = min(qq_expt)
  pMin, qMax = plotPQYieldSurfaceSim(uda_path,
                                     analytical_times,
                                     p_min_expt=pp_expt_min,
                                     q_min_expt=qq_expt_min)

  plt.title('Uniaxial strain SHPB ' + file_name)

  # Create output file name
  png_file_name = file_name + '_yield_surface'
  savePNG(os.path.join(save_path, png_file_name), '1280x960')
  #plt.show()

  #---------------------------------------------------------------------------------
  # Plot experimental and simulation data as a function of time
  analytical_times = np.array(analytical_times) * 1.0e6
  times = np.array(times) * 1.0e6
  time_expt = np.array(time_expt) * 1.0e6
  fig3 = plt.figure(3)
  plt.clf()
  plt.subplots_adjust(right=0.75)
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  plotExptDataSigmaTime(fig3, analytical_times, time_expt, sigma_a_expt,
                        sigma_r_expt)
  plotSimDataSigmaTime(fig3, analytical_times, times, sigma_a_sim, sigma_r_sim,
                       sigma_ar_sim)
  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  plt.xlabel(str_to_mathbf('Time (micro sec)'))
  plt.ylabel(str_to_mathbf('Axial and Radial Stress (MPa)'))
  plt.grid(True)
  plt.legend(loc='best', prop={'size': 10})
  plt.title('Uniaxial strain SHPB ' + file_name)
  png_file_name = file_name + '_sigma_time'
  savePNG(os.path.join(save_path, png_file_name), '1280x960')

  #---------------------------------------------------------------------------------
  # Plot experimental and simulation data as a function of time
  #analytical_times = np.array(analytical_times)*1.0e6
  #times = np.array(times)*1.0e6
  #time_expt = np.array(time_expt)*1.0e6
  pp_expt = np.array(pp_expt) * (-1.0)
  qq_expt = np.array(qq_expt) * (-1.0)
  pp_sim = np.array(pp_sim) * (-1.0)
  qq_sim = np.array(qq_sim) * (-1.0)
  fig4 = plt.figure(4)
  plt.clf()
  plt.subplots_adjust(right=0.75)
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  plotExptDataPQTime(fig4, analytical_times, time_expt, pp_expt, qq_expt)
  plotSimDataPQTime(fig4, analytical_times, times, pp_sim, qq_sim)
  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  plt.xlabel(str_to_mathbf('Time (micro sec)'))
  plt.ylabel(str_to_mathbf('Stress Invariants (MPa)'))
  plt.grid(True)
  plt.legend(loc='best', prop={'size': 10})
  plt.title('Uniaxial strain SHPB' + file_name)
  png_file_name = file_name + '_pq_time'
  savePNG(os.path.join(save_path, png_file_name), '1280x960')

  plt.show()


def test01_postProc(uda_path, save_path, **kwargs):
  print("Post Processing Test: 01 - Uniaxial Compression With Rotation")

  # Read the experimental data
  # file_name = 'Arena052212-035.expt'
  # file_name_alt = 'Arena052212-026.expt'
  # time_expt, sigma_a_expt, sigma_r_expt, pp_expt, qq_expt = readExptStressData(
  #     uda_path, file_name)
  # time_expt_alt, sigma_a_expt_alt, sigma_r_expt_alt, pp_expt_alt, qq_expt_alt = readExptStressData(
  #     uda_path, file_name_alt)

  # Read the simulation data
  times, sigmas, sigma_a_sim, sigma_r_sim, sigma_ar_sim, pp_sim, qq_sim = readSimStressData(
      uda_path)

  print(times)
  print(sigmas)

  # Get the model parameters
  material_dict = get_yield_surface(uda_path)
  param_text = material_dict['material string']

  # Set up time points
  #analytical_times = [0.0, 0.1, 0.2, 0.5, 0.7, 1.0]
  analytical_times = np.linspace(0.0, 0.001, 15)

  # Get snapshots of pq data (both expt and sim)
  # t_expt_snap, p_expt_snap = getDataTimeSnapshots(analytical_times, time_expt,
  #                                                 pp_expt)
  # t_expt_snap, q_expt_snap = getDataTimeSnapshots(analytical_times, time_expt,
  #                                                 qq_expt)
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
  plt.plot(times, Sxx, '-r', label='Uintah')
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
  plt.title('AreniscaTest 01:\nUniaxial Strain Compression')
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

  # Plot p vs. q simulation results
  eqShear_vs_meanStress(pp_sim, qq_sim)

  # Plot filled circles at time snapshots
  for ii in range(0, len(t_sim_snap) - 1):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(t_sim_snap))
    plt.plot(p_sim_snap[ii], q_sim_snap[ii], 'o', color=plt_color)

  # Plot the experimental data
  #line1 = plt.plot(pp_expt, qq_expt, '--b', linewidth=2, label='Expt. (35)')
  #line2 = plt.plot(pp_expt_alt,
  #                 qq_expt_alt,
  #                 '--m',
  #                 linewidth=2,
  #                 label='Expt. (26)')

  # Plot filled circles at time snapshots
  #for ii in range(0, len(t_expt_snap)):

    # Choose the Paired colormap
    #plt_color = cm.Paired(float(ii) / len(t_expt_snap))
    #plt.plot(p_expt_snap[ii], q_expt_snap[ii], 'v', color=plt_color)

  # Plot yield surfaces
  #pp_expt_min = min(pp_expt)
  #qq_expt_min = min(qq_expt)
  pMin, qMax = plotPQYieldSurfaceSim(uda_path,
                                     analytical_times,
                                     p_min_expt=0,
                                     q_min_expt=0)

  plt.title(
      'AreniscaTest: Dry Mason Sand: Uniaxial strain compression\n Yield surface evolution'
  )
  savePNG(save_path + '/Test01_yield_surface', '1280x960')
  #plt.show()

  #---------------------------------------------------------------------------------
  # Plot experimental and simulation data as a function of time
  fig3 = plt.figure(3)
  plt.clf()
  plt.subplots_adjust(right=0.75)
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  #plotExptDataSigmaTime(fig3, analytical_times, time_expt, sigma_a_expt,
  #                      sigma_r_expt)
  #plotExptDataSigmaTime(fig3, analytical_times, time_expt_alt, sigma_a_expt_alt, sigma_r_expt_alt)
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
  #plotExptDataPQTime(fig4, analytical_times, time_expt, pp_expt, qq_expt)
  #plotExptDataPQTime(fig4, analytical_times, time_expt_alt, pp_expt_alt, qq_expt_alt)
  plotSimDataPQTime(fig4, analytical_times, times, pp_sim, qq_sim)
  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  plt.xlabel(str_to_mathbf('Time (sec)'))
  plt.ylabel(str_to_mathbf('Stress (MPa)'))
  plt.grid(True)
  plt.legend(loc='best', prop={'size': 10})
  savePNG(save_path + '/Test01_pq_time', '1280x960')

  plt.show()


def test02_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 02 - Vertex Treatment")
  times, sigmas = get_pStress(uda_path)
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

  #Analytical Solutions
  #Drucker-Prager constants
  r0 = 50.0
  z0 = 50.0 * sqrtThree
  #Solution From Brannon Leelavanichkul paper
  analytical_times = [0.0, 1.0, threeHalves, 2.0, 5.0 / 2.0, 3.0]
  analytical_S11 = np.array([
      0, -850.0 / 3.0, (-50.0 / 3.0) * (9.0 + 4.0 * np.sqrt(6.0)),
      (-50.0 / 3.0) * (9.0 + 4.0 * np.sqrt(6.0)),
      (50.0 / 3.0) * (2.0 * np.sqrt(6) - 3.0),
      160.0 * np.sqrt(twoThirds) - 110.0
  ])
  analytical_S22 = np.array([
      0, -850.0 / 3.0, (50.0 / 3.0) * (2.0 * np.sqrt(6.0) - 9.0),
      (50.0 / 3.0) * (2.0 * np.sqrt(6.0) - 9.0),
      (-50.0 / 3.0) * (3.0 + np.sqrt(6.0)),
      (-10.0 / 3.0) * (33.0 + 8.0 * np.sqrt(6.0))
  ])
  analytical_S33 = np.array([
      0, -850.0 / 3.0, (50.0 / 3.0) * (2.0 * np.sqrt(6.0) - 9.0),
      (50.0 / 3.0) * (2.0 * np.sqrt(6.0) - 9.0),
      (-50.0 / 3.0) * (3.0 + np.sqrt(6.0)),
      (-10.0 / 3.0) * (33.0 + 8.0 * np.sqrt(6.0))
  ])
  analytical_mean = (analytical_S11 + analytical_S22 + analytical_S33) / 3.0
  analytical_I1 = analytical_S11 + analytical_S22 + analytical_S33
  tmp = (1.0 / 3.0) * analytical_I1
  analytical_s1 = analytical_S11 - tmp
  analytical_s2 = analytical_S22 - tmp
  analytical_s3 = analytical_S33 - tmp
  analytical_J2 = (1.0 / 2.0) * (pow(analytical_s1, 2) + pow(analytical_s2, 2) +
                                 pow(analytical_s3, 2))
  analytical_J3 = (1.0 / 3.0) * (pow(analytical_s1, 3) + pow(analytical_s2, 3) +
                                 pow(analytical_s3, 3))
  analytical_z = analytical_I1 / sqrtThree
  analytical_q = []
  for idx, J2 in enumerate(analytical_J2):
    J3 = analytical_J3[idx]
    analytical_q.append(sign(sqrtThree * np.sqrt(J2), J3))
  #Drucker-Prager yield surface
  yield_zs = np.array([z0, min(analytical_z)])
  yield_rs = r0 / z0 * (
      (get_yield_surface(uda_path)['PEAKI1'] / sqrtThree) - yield_zs)
  yield_ps = yield_zs * (sqrtThree / 3.0)
  yield_qs = yield_rs * np.sqrt(threeHalves)

  # Get the volumetric plastic strain
  times, plasticStrainVol = get_pPlasticStrainVol(uda_path)
  times, pCapX = get_capX(uda_path)
  times, pKappa = get_pKappa(uda_path)
  times, pZeta = get_zeta(uda_path)
  ev_p_list = []
  capX_list = []
  kappa_list = []
  zeta_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        ev_p_list.append(plasticStrainVol[ii])
        capX_list.append(pCapX[ii])
        #kappa_list.append(pKappa[ii])
        zeta_list.append(pZeta[ii])

  #print("ev_p = ", ev_p_list)
  #print("cap_X = ", capX_list)
  #print("kappa = ", kappa_list)
  #print("zeta = ", zeta_list)
  #print(sorted(set(ev_p_list)))
  ev_p_list_new = list(sorted(set(ev_p_list)))

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
  #plt.plot(analytical_mean,analytical_q,'-g',linewidth=lineWidth+1,label='Analytical')
  #plt.plot(yield_ps,yield_qs,'--k',linewidth=lineWidth+2,label='Yield surface')
  #plt.plot(yield_ps,-yield_qs,'--k',linewidth=lineWidth+2)
  #eqShear_vs_meanStress(ps,qs,(-300,60),(-300,300))
  #eqShear_vs_meanStress(ps,qs,(ps_min, ps_max),(qs_min, -qs_min))

  eqShear_vs_meanStress(ps, qs)

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('AreniscaTest 02:\nVertex Treatment (plot a)')
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
  plt.plot(analytical_times, analytical_S33, '-g', linewidth=lineWidth + 2)
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
  plt.plot(analytical_times,
           analytical_S11,
           '-g',
           linewidth=lineWidth + 2,
           label='Analytical')
  plt.plot(times, np.array(Sxx), '-r', label='Uintah')
  #Add Yield Surface
  #Add Analytical
  plt.legend()
  plt.setp(ax1.get_xticklabels(), visible=False)
  plt.ylabel(str_to_mathbf('\sigma_{xx} (MPa)'))
  plt.title('AreniscaTest 02:\nVertex Treatment (plot b)')
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  ax1.set_xlim(0, endT)
  #ax1.set_ylim(-400,100)
  #ax1.set_yticks([-400,-300,-200,-100,0,100])
  plt.grid(True)
  #Sigma yy
  ax2 = plt.subplot(312, sharex=ax3)
  plt.plot(analytical_times, analytical_S22, '-g', linewidth=lineWidth + 2)
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

  if SHOW_ON_MAKE:
    plt.show()


def test03_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 03 - Uniaxial Strain Without Hardening")
  times, sigmas = get_pStress(uda_path)
  ps_unscaled, qs_unscaled = get_ps_and_qs(sigmas)
  material_dict = get_yield_surface(uda_path)
  PEAKI1 = material_dict['PEAKI1']
  J2Yield = J2_at_Yield(uda_path)
  q_yield = np.sqrt(3.0 * J2Yield)

  # Scale the data
  ps = []
  qs = []
  for val in ps_unscaled:
    ps.append(val * 1.0e-6)

  for val in qs_unscaled:
    qs.append(val * 1.0e-6)

  PEAKI1 = PEAKI1 * 1.0e-6
  q_yield = q_yield * 1.0e-6
  #print('J2Yield : ',J2Yield)
  #print('q_yield : ',q_yield)

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
  # Get the volumetric plastic strain
  times, plasticStrainVol = get_pPlasticStrainVol(uda_path)
  times, pCapX = get_capX(uda_path)
  #times, pKappa = get_pKappa(uda_path)
  times, pZeta = get_zeta(uda_path)
  ev_p_list = []
  capX_list = []
  kappa_list = []
  zeta_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        ev_p_list.append(plasticStrainVol[ii])
        capX_list.append(pCapX[ii])
        #kappa_list.append(pKappa[ii])
        zeta_list.append(pZeta[ii])

  #print("ev_p = ", ev_p_list)
  #print("cap_X = ", capX_list)
  #print(kappa_list)
  #print("zeta = ", zeta_list)
  #print(sorted(set(ev_p_list)))
  ev_p_list_new = list(sorted(set(ev_p_list)))

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
  plt.title('AreniscaTest 03:\nUniaxial Strain Without Hardening')

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.plot(Xlims, (q_yield, q_yield),
           '--k',
           linewidth=lineWidth + 1,
           label='Initial yield surface')
  plt.plot(Xlims, (-q_yield, -q_yield), '--k', linewidth=lineWidth + 1)
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  #plt.legend()
  savePNG(save_path + '/Test03_verificationPlot', '1280x960')
  if SHOW_ON_MAKE:
    plt.show()


def test04_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 04 - Curved Yield Surface")
  times, sigmas = get_pStress(uda_path)
  ps_unscaled, qs_unscaled = get_ps_and_qs(sigmas)

  # Scale the data
  ps = []
  qs = []
  for val in ps_unscaled:
    ps.append(val * 1.0e-6)

  for val in qs_unscaled:
    qs.append(val * 1.0e-6)

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
  # Get the volumetric plastic strain
  times, plasticStrainVol = get_pPlasticStrainVol(uda_path)
  times, pCapX = get_capX(uda_path)
  #times, pKappa = get_pKappa(uda_path)
  times, pZeta = get_zeta(uda_path)
  ev_p_list = []
  capX_list = []
  kappa_list = []
  zeta_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        ev_p_list.append(plasticStrainVol[ii])
        capX_list.append(pCapX[ii])
        #kappa_list.append(pKappa[ii])
        zeta_list.append(pZeta[ii])

  #print("ev_p = ", ev_p_list)
  #print("cap_X = ", capX_list)
  #print(kappa_list)
  #print("zeta = ", zeta_list)
  #print(sorted(set(ev_p_list)))
  ev_p_list_new = list(sorted(set(ev_p_list)))

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

  plt.title('AreniscaTest 04:\nCurved Yield Surface')
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  #Add Analytical
  #plt.legend()
  savePNG(save_path + '/Test04_verificationPlot', '1280x960')
  if SHOW_ON_MAKE:
    plt.show()


def test05_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 05 - Hydrostatic Compression Fixed Cap")
  times, sigmas = get_pStress(uda_path)
  ps_unscaled, qs_unscaled = get_ps_and_qs(sigmas)

  # Scale the data
  ps = []
  qs = []
  for val in ps_unscaled:
    ps.append(val * 1.0e-6)

  for val in qs_unscaled:
    qs.append(val * 1.0e-6)

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
  # Get the volumetric plastic strain
  times, plasticStrainVol = get_pPlasticStrainVol(uda_path)
  times, pCapX = get_capX(uda_path)
  #times, pKappa = get_pKappa(uda_path)
  times, pZeta = get_zeta(uda_path)
  ev_p_list = []
  capX_list = []
  kappa_list = []
  zeta_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        ev_p_list.append(plasticStrainVol[ii])
        capX_list.append(pCapX[ii])
        #kappa_list.append(pKappa[ii])
        zeta_list.append(pZeta[ii])

  #print("ev_p = ", ev_p_list)
  #print("cap_X = ", capX_list)
  #print(kappa_list)
  #print("zeta = ", zeta_list)
  #print(sorted(set(ev_p_list)))
  ev_p_list_new = list(sorted(set(ev_p_list)))

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

  plt.title('AreniscaTest 05:\nHydrostatic Compression Fixed Cap')
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  #Add Analytical
  #plt.legend()
  savePNG(save_path + '/Test05_verificationPlot', '1280x960')
  if SHOW_ON_MAKE:
    plt.show()


def test06_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 06 - Uniaxial Strain Cap Evolution")
  times, sigmas = get_pStress(uda_path)
  ps_unscaled, qs_unscaled = get_ps_and_qs(sigmas)

  # Scale the data
  ps = []
  qs = []
  for val in ps_unscaled:
    ps.append(val * 1.0e-6)

  for val in qs_unscaled:
    qs.append(val * 1.0e-6)

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
  # Get the volumetric plastic strain
  times, plasticStrainVol = get_pPlasticStrainVol(uda_path)
  times, pCapX = get_capX(uda_path)
  #times, pKappa = get_pKappa(uda_path)
  times, pZeta = get_zeta(uda_path)
  ev_p_list = []
  capX_list = []
  kappa_list = []
  zeta_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        ev_p_list.append(plasticStrainVol[ii])
        capX_list.append(pCapX[ii])
        #kappa_list.append(pKappa[ii])
        zeta_list.append(pZeta[ii])

  #print("ev_p = ", ev_p_list)
  #print("cap_X = ", capX_list)
  #print(kappa_list)
  #print("zeta = ", zeta_list)
  #print(sorted(set(ev_p_list)))
  ev_p_list_new = list(sorted(set(ev_p_list)))

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
  plt.title('AreniscaTest 06:\nUniaxial Strain Cap Evolution')

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  #Add Analytical
  #plt.legend()
  savePNG(save_path + '/Test06_verificationPlot', '1280x960')
  if SHOW_ON_MAKE:
    plt.show()


def test07_postProc(uda_path, save_path, **kwargs):

  #Extract stress history
  print("Post Processing Test: 07 - Hydrostatic Compression Cap Evolution")
  times, sigmas = get_pStress(uda_path)
  I1s = []
  for sigma in sigmas:
    I1s.append(sigma_I1(sigma * 1.0e-6))

  # Find min/max values
  I1s_min = min(I1s)
  I1s_max = max(I1s)
  print("I1s_min = ", I1s_min)
  print("I1s_max = ", I1s_max)

  #analytical_times = [0.0,1.0]
  analytical_times = times

  # Get the internal variables
  ev_e, ev_p, capX, kappa, zeta, times = getAllInternalVariables(uda_path)

  # Get the material properties
  material_dict = get_yield_surface(uda_path)
  P3 = material_dict['P3']
  porosity = 1 - np.exp(-(P3 + np.array(ev_p)))

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  ##Plot a
  I1lims = (I1s_min, I1s_max)
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')

  ax1 = eqShear_vs_meanStress(I1s, porosity)
  plt.title('AreniscaTest 07:\nHydrostatic Compression with Fixed Cap')
  plt.ylabel(str_to_mathbf('Porosity'))
  plt.xlabel(str_to_mathbf('I_{1}:first invariant of stress tensor (MPa)'))
  #plt.show()

  plot_crush_curve(uda_path, I1lims)
  #ax1.set_xticks([-9000,-7000,-5000,-3000,-1000,0])
  #ax1.set_xticks([-8000,-6000,-4000,-2000,0])
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  plt.legend()
  savePNG(save_path + '/Test07_verificationPlot', '1280x960')
  plt.show()
  if SHOW_ON_MAKE:
    plt.show()


def test08_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 08 - Loading/Unloading")
  times, sigmas = get_pStress(uda_path)
  I1s = []
  ps = []
  for sigma in sigmas:
    I1s.append(sigma_I1(sigma * 1.0e-6))
    ps.append(sigma_I1(sigma) / 3.0 * 1.0e-6)

  # Find min/max values
  I1s_min = min(I1s)
  I1s_max = max(I1s)
  ps_min = min(ps)
  ps_max = max(ps)
  print("I1s_min = ", I1s_min)
  print("I1s_max = ", I1s_max)
  print("ps_min = ", ps_min)
  print("ps_max = ", ps_max)

  analytical_times = [0.0, 1.0, 2.0, 3.0, 4.0]
  # Get the volumetric plastic strain
  times, plasticStrainVol = get_pPlasticStrainVol(uda_path)
  times, elasticStrainVol = get_pElasticStrainVol(uda_path)
  times, pCapX = get_capX(uda_path)
  times, pZeta = get_zeta(uda_path)
  ev_p_list = []
  capX_list = []
  kappa_list = []
  zeta_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 1.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        ev_p_list.append(plasticStrainVol[ii])
        capX_list.append(pCapX[ii])
        #kappa_list.append(pKappa[ii])
        zeta_list.append(pZeta[ii])

  #print("ev_p = ", ev_p_list)
  #print("cap_X = ", capX_list)
  #print(kappa_list)
  #print("zeta = ", zeta_list)
  #print(sorted(set(ev_p_list)))
  ev_p_list_new = list(sorted(set(ev_p_list)))

  totalStrainVol = np.array(elasticStrainVol) + np.array(plasticStrainVol)
  material_dict = get_yield_surface(uda_path)
  P3 = material_dict['P3']
  porosity = 1 - np.exp(-(P3 + np.array(plasticStrainVol)))

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
  print(len(times))
  print(len(ps))
  ax1 = eqShear_vs_meanStress(times, -np.array(ps))

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('AreniscaTest 08:\nLoading/Unloading (plot a)')
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

  plt.title('AreniscaTest 08:\nLoading/Unloading (plot b)')
  plt.ylabel(str_to_mathbf('Total Volumetric Strain, \epsilon_{v}'))
  plt.xlabel(str_to_mathbf('Time (s)'))
  ax2.xaxis.set_major_formatter(int_formatter)
  ax2.yaxis.set_major_formatter(int_formatter)
  ax2.tick_params(axis='both', labelsize='small')
  savePNG(save_path + '/Test08_verificationPlot_b', '1280x960')

  ##Plot c
  I1lims = (I1s_min, I1s_max)
  plt.figure(3)
  plt.clf()
  ax3 = plt.subplot(111)
  plt.subplots_adjust(right=0.75, left=0.15)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  eqShear_vs_meanStress(I1s, porosity, I1lims, (0, 1.25))
  plt.title('AreniscaTest 08:\nLoading/Unloading (plot c)')
  plt.ylabel(str_to_mathbf('Porosity'))
  plt.xlabel(str_to_mathbf('I_{1}:first invariant of stress tensor (Pa)'))
  plot_crush_curve(uda_path, I1lims)
  #ax1.set_xticks([-9000,-7000,-5000,-3000,-1000,0])
  #ax3.set_xticks([-10000,-7500,-5000,-2500,0,1000])
  ax3.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
  ax3.xaxis.set_major_formatter(exp_formatter)
  ax3.yaxis.set_major_formatter(int_formatter)
  ax3.tick_params(axis='both', labelsize='small')
  plt.legend()
  savePNG(save_path + '/Test08_verificationPlot_c', '1280x960')

  if SHOW_ON_MAKE:
    plt.show()


def test09_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 09 - Fluid Filled Pore Space")
  times, sigmas = get_pStress(uda_path)
  I1s = []
  ps = []
  for sigma in sigmas:
    I1s.append(sigma_I1(sigma) * 1.0e-6)
    ps.append(sigma_I1(sigma) / 3.0 * 1.0e-6)

  # Find min/max values
  I1s_min = min(I1s)
  I1s_max = max(I1s)
  ps_min = min(ps)
  ps_max = max(ps)
  print("I1s_min = ", I1s_min)
  print("I1s_max = ", I1s_max)
  print("ps_min = ", ps_min)
  print("ps_max = ", ps_max)

  analytical_times = [0.0, 1.0, 2.0, 3.0, 4.0]
  # Get the volumetric plastic strain
  times, plasticStrainVol = get_pPlasticStrainVol(uda_path)
  times, elasticStrainVol = get_pElasticStrainVol(uda_path)
  times, pCapX = get_capX(uda_path)
  times, pZeta = get_zeta(uda_path)
  ev_p_list = []
  capX_list = []
  kappa_list = []
  zeta_list = []
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times) - 1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 1.0e-4 and jj > 0):
        #print("ii = " , ii, " tt = " , tt, "jj = " , jj, " ta = " , ta )
        ev_p_list.append(plasticStrainVol[ii])
        capX_list.append(pCapX[ii])
        #kappa_list.append(pKappa[ii])
        zeta_list.append(pZeta[ii])

  #print("ev_p = ", ev_p_list)
  #print("cap_X = ", capX_list)
  #print(kappa_list)
  #print("zeta = ", zeta_list)
  #print(sorted(set(ev_p_list)))
  ev_p_list_new = list(sorted(set(ev_p_list)))

  totalStrainVol = np.array(elasticStrainVol) + np.array(plasticStrainVol)
  material_dict = get_yield_surface(uda_path)
  P3 = material_dict['P3']
  porosity = 1 - np.exp(-(P3 + np.array(plasticStrainVol)))

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
  ax1 = eqShear_vs_meanStress(times, -np.array(ps))
  plt.title('AreniscaTest 09:\nFluid EFfects (plot a)')
  plt.ylabel(str_to_mathbf('Pressure (Pa)'))
  plt.xlabel(str_to_mathbf('Time (s)'))
  ax1.xaxis.set_major_formatter(int_formatter)
  ax1.yaxis.set_major_formatter(exp_formatter)
  ax1.tick_params(axis='both', labelsize='small')
  savePNG(save_path + '/Test09_verificationPlot_a', '1280x960')

  ##Plot b
  plt.figure(2)
  plt.clf()
  ax2 = plt.subplot(111)
  plt.subplots_adjust(right=0.75, left=0.15)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  ax1 = eqShear_vs_meanStress(times, totalStrainVol)
  plt.title('AreniscaTest 09:\nFluid EFfects (plot b)')
  plt.ylabel(str_to_mathbf('Total Volumetric Strain, \epsilon_{v}'))
  plt.xlabel(str_to_mathbf('Time (s)'))
  ax2.xaxis.set_major_formatter(int_formatter)
  ax2.yaxis.set_major_formatter(int_formatter)
  ax2.tick_params(axis='both', labelsize='small')
  savePNG(save_path + '/Test09_verificationPlot_b', '1280x960')

  ##Plot c
  I1lims = (I1s_min, I1s_max)
  plt.figure(3)
  plt.clf()
  ax3 = plt.subplot(111)
  plt.subplots_adjust(right=0.75, left=0.15)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  eqShear_vs_meanStress(I1s, porosity, I1lims, (0, 1.25))
  plt.title('AreniscaTest 09:\nFluid EFfects (plot c)')
  plt.ylabel(str_to_mathbf('Porosity'))
  plt.xlabel(str_to_mathbf('I_{1}:first invariant of stress tensor (Pa)'))
  plot_crush_curve(uda_path, I1lims)
  #ax1.set_xticks([-9000,-7000,-5000,-3000,-1000,0])
  ax3.set_xticks([-10000, -7500, -5000, -2500, 0, 1000])
  ax3.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
  ax3.xaxis.set_major_formatter(exp_formatter)
  ax3.yaxis.set_major_formatter(int_formatter)
  ax3.tick_params(axis='both', labelsize='small')
  plt.legend()
  savePNG(save_path + '/Test09_verificationPlot_c', '1280x960')

  if SHOW_ON_MAKE:
    plt.show()


def test10_postProc(uda_path, save_path, **kwargs):

  BIG_FIGURE = True
  if 'WORKING_PATH' in kwargs:
    working_dir = kwargs['WORKING_PATH']
  else:
    print('\nERROR: need working directory to post process this problem')
    return

  #Extract stress history
  print(
      "Post Processing Test: 10 - Transient Stress Eigenvalues with Constant Eigenvectors"
  )
  times, sigmas = get_pStress(uda_path)
  Sxx = []
  Syy = []
  Szz = []
  for sigma in sigmas:
    Sxx.append(sigma[0][0] * 1.0e-6)
    Syy.append(sigma[1][1] * 1.0e-6)
    Szz.append(sigma[2][2] * 1.0e-6)

  # Find min/max values
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

  #Analytical solution
  material_dict = get_yield_surface(uda_path)
  def_times, Fs = get_defTable(uda_path, working_dir)
  tau_yield = material_dict['PEAKI1']
  bulk_mod = material_dict['B0']
  shear_mod = material_dict['G0']

  analytical_times, analytical_sigmas, epsils = defTable_to_J2Solution(
      def_times, Fs, bulk_mod, shear_mod, tau_yield, num_substeps=10)

  analytical_Sxx = []
  analytical_Syy = []
  analytical_Szz = []
  for sigma in analytical_sigmas:
    analytical_Sxx.append(sigma[0][0] * 1.0e-6)
    analytical_Syy.append(sigma[1][1] * 1.0e-6)
    analytical_Szz.append(sigma[2][2] * 1.0e-6)

  ###PLOTTING
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  if BIG_FIGURE:
    plt.subplots_adjust(right=0.75)
    param_text = material_dict['material string']
    plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  else:
    plt.subplots_adjust(left=0.15, top=0.96, bottom=0.15, right=0.96)

  #analytical solution
  plt.plot(analytical_times,
           np.array(analytical_Sxx),
           ':r',
           linewidth=lineWidth + 2,
           label=str_to_mathbf('Analytical \sigma_{xx}'))
  plt.plot(analytical_times,
           np.array(analytical_Syy),
           '--g',
           linewidth=lineWidth + 2,
           label=str_to_mathbf('Analytical \sigma_{yy}'))
  plt.plot(analytical_times,
           np.array(analytical_Szz),
           '-.b',
           linewidth=lineWidth + 2,
           label=str_to_mathbf('Analytical \sigma_{zz}'))

  #simulation results
  plt.plot(times,
           np.array(Sxx),
           '-r',
           label=str_to_mathbf('Uintah \sigma_{xx}'))
  plt.plot(times,
           np.array(Syy),
           '-g',
           label=str_to_mathbf('Uintah \sigma_{yy}'))
  plt.plot(times,
           np.array(Szz),
           '-b',
           label=str_to_mathbf('Uintah \sigma_{zz}'))

  ax1.set_xlim(0, 2.25)
  formatter_int = ticker.FormatStrFormatter('$\mathbf{%g}$')
  ax1.xaxis.set_major_formatter(formatter_int)
  ax1.yaxis.set_major_formatter(formatter_int)
  #labels
  plt.grid(True)
  plt.xlabel(str_to_mathbf('Time (s)'))
  plt.ylabel(str_to_mathbf('Stress (MPa)'))
  if BIG_FIGURE:
    plt.legend(loc='upper right', bbox_to_anchor=(1.38, 1.12))
    plt.title(
        'AreniscaTest 10:\nTransient Stress Eigenvalues with Constant Eigenvectors'
    )
    savePNG(save_path + '/Test10_verificationPlot', '1280x960')
  else:
    tmp = plt.rcParams['legend.fontsize']
    plt.rcParams['legend.fontsize'] = 'x-small'
    plt.legend(loc=7)
    savePNG(save_path + '/Test10_verificationPlot', '640x480')
    plt.rcParams['legend.fontsize'] = tmp

  if SHOW_ON_MAKE:
    plt.show()


def test11_postProc(uda_path, save_path, **kwargs):
  if 'WORKING_PATH' in kwargs:
    working_dir = kwargs['WORKING_PATH']
  else:
    print('\nERROR: need working directory to post process this problem')

  #Extract stress and strain history
  print("Post Processing Test: 11 - Uniaxial Strain J2 Plasticity")
  times, sigmas = get_pStress(uda_path)
  times, epsils = get_epsilons(uda_path)
  exx = []
  eyy = []
  ezz = []
  for epsil in epsils:
    exx.append(epsil[0][0])
    eyy.append(epsil[1][1])
    ezz.append(epsil[2][2])
  Sxx = []
  Syy = []
  Szz = []
  for sigma in sigmas:
    Sxx.append(sigma[0][0] * 1.0e-6)
    Syy.append(sigma[1][1] * 1.0e-6)
    Szz.append(sigma[2][2] * 1.0e-6)

  #Analytical solution
  material_dict = get_yield_surface(uda_path)
  def_times, Fs = get_defTable(uda_path, working_dir)
  tau_yield = material_dict['PEAKI1'] * material_dict['FSLOPE']
  #tau_yield = material_dict['PEAKI1']
  bulk_mod = material_dict['B0']
  shear_mod = material_dict['G0']

  analytical_times, analytical_sigmas, epsils = defTable_to_J2Solution(
      def_times, Fs, bulk_mod, shear_mod, tau_yield, num_substeps=1000)

  analytical_e11 = []
  analytical_e22 = []
  analytical_e33 = []
  for epsil in epsils:
    analytical_e11.append(epsil[0][0])
    analytical_e22.append(epsil[1][1])
    analytical_e33.append(epsil[2][2])

  analytical_Sxx = []
  analytical_Syy = []
  analytical_Szz = []
  for sigma in analytical_sigmas:
    analytical_Sxx.append(sigma[0][0] * 1.0e-6)
    analytical_Syy.append(sigma[1][1] * 1.0e-6)
    analytical_Szz.append(sigma[2][2] * 1.0e-6)

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  plt.title('AreniscaTest 11:\nUniaxial Strain J2 Plasticity (plot a)')
  plt.plot(np.array(analytical_e11),
           np.array(analytical_Sxx),
           '--g',
           linewidth=lineWidth + 1,
           label=str_to_mathbf('Analytical'))
  plt.plot(np.array(exx), np.array(Sxx), '-r', label=str_to_mathbf('Uintah'))
  plt.xlabel(str_to_mathbf('\epsilon_{A}'))
  plt.ylabel(str_to_mathbf('\sigma_{A} (MPa)'))
  plt.legend()
  savePNG(save_path + '/Test11_verificationPlot_a', '1280x960')

  plt.figure(2)
  plt.clf()
  ax2 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  ax2.xaxis.set_major_formatter(formatter)
  ax2.yaxis.set_major_formatter(formatter)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  plt.title('AreniscaTest 11:\nUniaxial Strain J2 Plasticity (plot b)')
  plt.plot(np.array(analytical_e11),
           np.array(analytical_Syy),
           '--g',
           linewidth=lineWidth + 1,
           label=str_to_mathbf('Analytical'))
  plt.plot(np.array(exx), np.array(Syy), '-r', label=str_to_mathbf('Uintah'))
  plt.xlabel(str_to_mathbf('\epsilon_{A}'))
  plt.ylabel(str_to_mathbf('\sigma_{L} (MPa)'))
  plt.legend()
  savePNG(save_path + '/Test11_verificationPlot_b', '1280x960')

  plt.figure(3)
  plt.clf()
  ax3 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  ax3.xaxis.set_major_formatter(formatter)
  ax3.yaxis.set_major_formatter(formatter)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  plt.title('AreniscaTest 11:\nUniaxial Strain J2 Plasticity (plot c)')
  plt.plot(analytical_times,
           np.array(analytical_e11),
           '-g',
           linewidth=lineWidth + 1,
           label=str_to_mathbf('Analytical \epsilon_{xx}'))
  plt.plot(analytical_times,
           np.array(analytical_e22),
           '-r',
           linewidth=lineWidth + 1,
           label=str_to_mathbf('Analytical \epsilon_{yy}'))
  plt.plot(analytical_times,
           np.array(analytical_e33),
           '-b',
           linewidth=lineWidth + 1,
           label=str_to_mathbf('Analytical \epsilon_{zz}'))
  plt.legend()
  plt.xlabel(str_to_mathbf('Time (s)'))
  plt.ylabel(str_to_mathbf('\epsilon'))
  savePNG(save_path + '/Test11_verificationPlot_c', '1280x960')

  plt.figure(4)
  plt.clf()
  ax4 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  ax4.xaxis.set_major_formatter(formatter)
  ax4.yaxis.set_major_formatter(formatter)
  param_text = material_dict['material string']
  plt.figtext(0.77, 0.70, param_text, ha='left', va='top', size='x-small')
  plt.title('AreniscaTest 11:\nUniaxial Strain J2 Plasticity (plot d)')
  plt.plot(analytical_times,
           np.array(analytical_Sxx),
           '-g',
           linewidth=lineWidth + 1,
           label=str_to_mathbf('Analytical \sigma_{xx}'))
  plt.plot(analytical_times,
           np.array(analytical_Syy),
           '-r',
           linewidth=lineWidth + 1,
           label=str_to_mathbf('Analytical \sigma_{yy}'))
  plt.plot(analytical_times,
           np.array(analytical_Szz),
           '-b',
           linewidth=lineWidth + 1,
           label=str_to_mathbf('Analytical \sigma_{zz}'))
  plt.legend()
  plt.xlabel(str_to_mathbf('Time (s)'))
  plt.ylabel(str_to_mathbf('\sigma (MPa)'))
  savePNG(save_path + '/Test11_verificationPlot_d', '1280x960')

  if SHOW_ON_MAKE:
    plt.show()


def test12_postProc(uda_path, save_path, **kwargs):
  #Extract stress history
  print("Post Processing Test: 12 - Nonlinear Elasticity")
  times, sigmas = get_pStress(uda_path)
  pressure = []
  for sigma in sigmas:
    pressure.append(-sigma_I1(sigma) / 3.0 * 1.0e-6)
  times, plasticStrainVol = get_pPlasticStrainVol(uda_path)
  times, elasticStrainVol = get_pElasticStrainVol(uda_path)
  totalStrainVol = -np.array(elasticStrainVol) - np.array(plasticStrainVol)

  # Find min/max values
  ev_min = min(totalStrainVol)
  ev_max = max(totalStrainVol)
  print("ev_min = ", ev_min)
  print("ev_max = ", ev_max)

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  ##Plot a
  evlims = (ev_min, ev_max)
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  #param_text = material_dict['material string']
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')
  #ax1=eqShear_vs_meanStress(I1s,porosity,I1lims,(0,0.6))
  plt.plot(totalStrainVol, pressure, '-b', label='Arenisca')
  plt.title('AreniscaTest 12:\nNonlinear Elasticity')
  plt.ylabel(str_to_mathbf('p: pressure (MPa)'))
  plt.xlabel(str_to_mathbf('ev: compressive volumetric strain'))
  #ax1.set_xticks([0,0.005,0.010,0.015,0.020,0.025])
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  plt.legend()
  savePNG(save_path + '/Test12_verificationPlot', '1280x960')
  if SHOW_ON_MAKE:
    plt.show()


def test13_postProc(uda_path, save_path, **kwargs):
  COLORS = ['Black', 'Blue', 'Magenta', 'Red', 'Green']
  if 'WORKING_PATH' in kwargs:
    working_dir = kwargs['WORKING_PATH']
  else:
    print('\nERROR: need working directory to post process this problem')
    return

  #Plot Constants
  Xlims = (-450, 50)
  Ylims = (-100, 100)
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  plt.figure(1)
  plt.hold(True)
  plt.clf()

  material_dict = get_yield_surface(uda_path)
  PEAKI1 = material_dict['PEAKI1']
  FSLOPE = material_dict['FSLOPE']
  STREN = material_dict['STREN']
  T1 = material_dict['T1']
  T2 = material_dict['T2']

  def_times, Fs = get_defTable(uda_path, working_dir)
  A = Fs[1][0][0]
  #As = Fs[10][0][0]
  K = material_dict['B0']
  G = material_dict['G0']
  C = K + (4.0 / 3.0) * G
  Y = STREN * 1.732
  YS = STREN

  #uniaxial strain (unscaled)
  analytical_exx = [
      0.0,
      (Y / (2.0 * G)),
      np.log(A),
  ]
  analytical_Sxx = [
      0.0,
      (C * Y) / (2.0 * G) * 1.0e-6,
      (((C - K) * Y) / (2 * G) + K * np.log(A)) * 1.0e-6,
  ]

  #uniaxial strain (scaled)
  #analytical_exx = np.array([0.0,
  #(Y/(2.0*G)),
  #np.log(A),
  #np.log(A)-(Y)/(G),
  #0.0
  #])/(Y/(2.0*G))

  #analytical_Sxx = np.array([0.0,
  #(C*Y)/(2.0*G),
  #((C-K)*Y)/(2*G)+K*np.log(A),
  #K*np.log(A)-((C+K)*Y)/(2*G),
  #(K-C)*Y/(2*G)
  #])/((C*Y)/(2.0*G))

  #pure shear (unscaled)
  #analytical_exx = np.array([0.0,
  #          (YS/(2.0*G)),
  #          np.log(As),
  #          ])
  #analytical_Sxx = np.array([0.0,
  #          (YS),
  #          (YS),
  #          ])

  #Extract stress history
  print("Post Processing Test: 13 ")
  times, sigmas = get_pStress(uda_path)
  times, epsils = get_epsilons(uda_path)
  exx = []
  eyy = []
  ezz = []
  exy = []
  for epsil in epsils:
    exx.append(epsil[0][0])
    eyy.append(epsil[1][1])
    ezz.append(epsil[2][2])
    exy.append(epsil[0][1])
  Sxx = []
  Syy = []
  Szz = []
  Sxy = []
  for sigma in sigmas:
    Sxx.append(sigma[0][0] * 1.0e-6)
    Syy.append(sigma[1][1] * 1.0e-6)
    Szz.append(sigma[2][2] * 1.0e-6)
    Sxy.append(sigma[0][1] * 1.0e-6)

  scaled_exx = ((2.0 * G) / Y) * np.array(exx)
  scaled_Sxx = ((2.0 * G) / (C * Y)) * np.array(Sxx) * 1.0e6
  scaled_Syy = ((2.0 * G) / (C * Y)) * np.array(Syy) * 1.0e6
  #S = np.array(Sxx) - np.array(Syy)
  S = np.array(Sxx)
  #E = np.array(exy)

  ###PLOTTING
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  #param_text = material_dict['material string']
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')
  eqShear_vs_meanStress(exx,
                        S,
                        LINE_LABEL='T1=' + format(T1, '1.3e') + ' T2=' +
                        format(T2, '1.3e'))
  #eqShear_vs_meanStress(E,S,LINE_LABEL = 'T1='+format(T1,'1.3e')+' T2='+format(T2,'1.3e'),COLOR=COLORS[idx])
  plt.plot(analytical_exx,
           analytical_Sxx,
           '--',
           color='Red',
           label='Analytical solution for rate independent case.')
  plt.title('AreniscaTest 13:')
  plt.ylabel(str_to_mathbf('\sigma_{xx}'))
  plt.xlabel(str_to_mathbf('\epsilon_{xx}'))
  #plt.ylabel(str_to_mathbf('\sigma_{xy}'))
  #plt.xlabel(str_to_mathbf('\epsilon_{xy}'))
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)
  plt.legend()
  savePNG(save_path + '/Test13_verificationPlot', '1280x960')
  if SHOW_ON_MAKE:
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
  times, ev_e_list = get_pElasticStrainVol(uda_path)
  times, ev_p_list = get_pPlasticStrainVol(uda_path)
  times, capX_list = get_capX(uda_path)
  times, kappa_list = get_pKappa(uda_path)
  time_list, zeta_list = get_zeta(uda_path)

  return ev_e_list, ev_p_list, capX_list, kappa_list, zeta_list, time_list


#---------------------------------------------------------------------------------
# Read the internal state variables from the UDA file:
#    1) volumetric plastic strain
#    2) capX
#    3) kappa
#    3) zeta
#---------------------------------------------------------------------------------
def getInternalVariables(uda_path, analytical_times):

  # Get the internal variables
  ev_e, ev_p, capX, kappa, zeta, times = getAllInternalVariables(uda_path)

  # Create a clean list containing the data as functions of time
  ev_e_list = []
  ev_p_list = []
  capX_list = []
  kappa_list = []
  zeta_list = []
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
        ev_e_list.append((1 - ss) * ev_e[ii - 1] + ss * ev_e[ii])
        ev_p_list.append((1 - ss) * ev_p[ii - 1] + ss * ev_p[ii])
        capX_list.append((1 - ss) * capX[ii - 1] + ss * capX[ii])
        kappa_list.append((1 - ss) * kappa[ii - 1] + ss * kappa[ii])
        zeta_list.append((1 - ss) * zeta[ii - 1] + ss * zeta[ii])
    #    break
    #else:
    #  continue

  #print("time = ", time_list)
  #print("ev_e = ", ev_e_list)
  #print("ev_p = ", ev_p_list)
  #print("cap_X = ", capX_list)
  #print("kappa = ", kappa_list)
  #print("zeta = ", zeta_list)

  return ev_e_list, ev_p_list, capX_list, kappa_list, zeta_list, time_list


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
  ev_e_list, ev_p_list, capX_list, kappa_list, zeta_list, time_list = getInternalVariables(
      uda_path, time_points)

  # Get material parameters
  material_dict = get_yield_surface(uda_path)
  PEAKI1 = material_dict['PEAKI1']
  FSLOPE = material_dict['FSLOPE']
  STREN = material_dict['STREN']
  YSLOPE = material_dict['YSLOPE']
  Pf0 = material_dict['P_f0']

  # Set up constants
  a1 = STREN
  a2 = (FSLOPE - YSLOPE) / (STREN - YSLOPE * PEAKI1)
  a3 = (STREN - YSLOPE * PEAKI1) * math.exp(-a2 * PEAKI1)
  a4 = YSLOPE

  # Find min capX
  minCapX = min(capX_list)
  #print("min(CapX) = ", minCapX)

  # Find pmin and qmin
  pmin = 1.0e8
  qmax = -1.0e8

  # Loop through time snapshots
  for ii, ev_p in enumerate(ev_p_list):

    # Get the internal variables
    capX = capX_list[ii]
    kappa = kappa_list[ii]
    zeta = zeta_list[ii]

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(ev_p_list))

    # Create an array of I1 values
    num_points = 100
    I1s = np.linspace(0.999 * minCapX - 3.0 * Pf0, PEAKI1 - 3.0 * Pf0,
                      num_points)
    J2s = []

    #J2 versus I1
    for I1 in I1s:

      # Add the fluid pressure
      I1f = I1 + 3.0 * Pf0

      # Kinematic hardening shift
      I1_minus_zeta = I1f - zeta

      # If I1 < capX or I1 > PEAKI1 move I1 to yield surface
      if (I1_minus_zeta < capX):
        I1_minus_zeta = 0.999999 * capX

      if (I1_minus_zeta > PEAKI1):
        I1_minus_zeta = 0.999999 * PEAKI1

      # Compute F_f
      Ff = a1 - a3 * math.exp(a2 * I1_minus_zeta) - a4 * (I1_minus_zeta)
      Ff_sq = Ff**2

      # Compute Fc
      Fc_sq = 1.0
      if (I1_minus_zeta < kappa):
        #print("Computing cap")
        ratio = (kappa - I1_minus_zeta) / (kappa - capX)
        Fc_sq = 1.0 - ratio**2

      # Compute J2
      J2 = Ff_sq * Fc_sq
      J2s.append(J2)

    xs = np.array(I1s) / 3.0 * 1.0e-6
    ys = np.sqrt(3.0 * np.array(J2s)) * 1.0e-6

    if (min(xs) < pmin):
      pmin = min(xs)

    if (max(ys) > qmax):
      qmax = max(ys)

    #print(xs.shape, ys.shape)
    #print("ii = ", ii)
    ev_e_str = "%.2g" % (-ev_e_list[ii])
    ev_p_str = "%.2g" % (-ev_p)
    time_str = "%.2g" % time_list[ii]
    #print("ev_p " , ev_p_str)
    #print("time " , time_str)

    #label_str = '$\epsilon_v^e$ = ' + ev_e_str + ' $\epsilon_v^p$ = ' + ev_p_str + ' $t$ = ' + time_str
    label_str = ' $\epsilon_v^p$ = ' + ev_p_str
    line1 = plt.plot(-xs, ys, '--b', linewidth=1, label=label_str)
    line2 = plt.plot(-xs, -ys, '--b', linewidth=1)
    plt.setp(line1, color=plt_color)
    plt.setp(line2, color=plt_color)
    plt.legend(loc=2, prop={'size': 8})

  print("pmin = ", pmin, " qmax = ", qmax)
  print("pmin_expt = ", p_min_expt, " qmin_expt = ", q_min_expt)
  if (pmin > p_min_expt):
    pmin = p_min_expt
  if (qmax < -q_min_expt):
    qmax = -q_min_expt
  axes = plt.gca()
  #axes.set_xlim([1.3*pmin, 1])
  #axes.set_ylim([-1.3*qmax, 1.3*qmax])
  axes.set_xlim([-1, -1.3 * pmin])
  axes.set_ylim([-1.3 * qmax, 1.3 * qmax])
  return pmin, qmax


def plotCrushCurveExptSim(uda_path, I1lims=[-10000, 0]):
  nPoints = 500
  material_dict = get_yield_surface(uda_path)
  P0 = material_dict['P0']
  P1 = material_dict['P1']
  P3 = material_dict['P3']

  # Analytical solution for porosity vs. X for piece-wise crush curve
  # compression
  I1sC = np.linspace(I1lims[0] * 1.0e6, P0, nPoints)
  porosityC = 1 - np.exp(-P3 * np.exp(P1 * (I1sC - P0)))
  # tension
  I1sT = np.linspace(P0, I1lims[1] * 1.0e6, nPoints)
  porosityT = 1 - np.exp(-(I1sT / P0)**(P0 * P1 * P3) - P3 + 1)

  # Scale down again
  I1_c = []
  for ii in I1sC:
    I1_c.append(ii * 1.0e-6)
  I1_t = []
  for ii in I1sT:
    I1_t.append(ii * 1.0e-6)

  plt.plot(I1_c,
           porosityC,
           '--g',
           linewidth=lineWidth + 1,
           label='Analytical crush curve - Compression')
  plt.hold(True)
  plt.plot(I1_t,
           porosityT,
           '--b',
           linewidth=lineWidth + 1,
           label='Analytical crush curve - Tension')
