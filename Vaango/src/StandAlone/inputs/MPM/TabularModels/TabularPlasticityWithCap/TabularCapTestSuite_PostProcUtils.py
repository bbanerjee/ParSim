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

from TabularCapYieldSurfaceUtils import *

SHOW_ON_MAKE = False

# Get the partextact executable path from the environment
partextract_exe = os.path.abspath(os.environ['PARTEXTRACT_EXE'])

#Useful constants
sqrtThree = np.sqrt(3.0)
twoThirds = 2.0/3.0
threeHalves = 3.0/2.0

#Set matplotlib defaults to desired values
#Set the legend to best fit
fontSize = 16

markers = None
plt.rcParams['legend.loc']='best'
#Set font size
plt.rcParams['mathtext.it'] = 'serif:bold'
plt.rcParams['mathtext.rm'] = 'serif:bold'
plt.rcParams['mathtext.sf'] = 'serif:bold'
plt.rcParams['font.size']=fontSize
plt.rcParams['font.weight']='bold'
plt.rcParams['axes.labelsize']='medium'
#plt.rcParams['axes.labelweight']='bold'
plt.rcParams['legend.fontsize']='medium'
#Set linewidth
lineWidth = 2
plt.rcParams['lines.linewidth']=lineWidth
#Set markersize
plt.rcParams['lines.markersize'] = 8
#Set padding for tick labels and size
plt.rcParams['xtick.major.pad']  = 12
plt.rcParams['ytick.major.pad']  = 8
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.size'] = 3

#resolution
plt.rcParams['figure.dpi']=120

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : fontSize}
rc('font', **font) 
rc('text', usetex=True)
 

def savePNG(name,size='1920x1080'):
  res = float(plt.rcParams['figure.dpi'])
  print(res)
  #Add Check for file already existing as name.png
  if size == '640x480':
    size = [640/res,480/res]
  if size == '1080x768':
    size = [1080/res,768/res]
  if size == '1152x768':
    size = [1152/res,768/res]
  if size == '1280x854':
    size = [1280/res,854/res]    
  if size == '1280x960':
    size = [1280/res,960/res]
  if size == '1920x1080':
    size = [1920/res,1080/res]
  #set the figure size for saving
  #plt.gcf().set_size_inches(size[0],size[1])
  #save at speciified resolution
  #plt.savefig(name+'.png', bbox_inches=0, dpi=plt.rcParams['figure.dpi']) 

  #size='960x960'
  print(size)
  plt.gcf().set_size_inches(size[0],size[1])
  plt.savefig(name+'.pdf', bbox_inches=0, dpi=plt.rcParams['figure.dpi']) 

def savePDF(name,size='1920x1080'):
  res = float(plt.rcParams['figure.dpi'])
  #Add Check for file already existing as name.png
  if size == '640x480':
    size = [640/res,480/res]
  if size == '1080x768':
    size = [1080/res,768/res]
  if size == '1152x768':
    size = [1152/res,768/res]
  if size == '1280x854':
    size = [1280/res,854/res]    
  if size == '1280x960':
    size = [1280/res,960/res]
  if size == '1920x1080':
    size = [1920/res,1080/res]
  #set the figure size for saving
  plt.gcf().set_size_inches(size[0],size[1])
  #save at speciified resolution
  plt.savefig(name+'.pdf', bbox_inches=0, dpi=plt.rcParams['figure.dpi']) 

def sign(x,y):
  if y>=0:
    return abs(x)
  else:
    return -abs(x)

def sigma_iso(sigma):
  return (np.trace(sigma)/3.0)*np.eye(3)
  
def sigma_dev(sigma):
  return sigma-sigma_iso(sigma)

def sigma_I1(sigma):
  #print(sigma)
  return sigma.trace()
  
def sigma_J2(sigma):
  return 0.5*np.dot(sigma_dev(sigma),sigma_dev(sigma)).trace()

def sigma_J3(sigma):
  return (1/3.0)*np.dot(np.dot(sigma_dev(sigma),sigma_dev(sigma)),sigma_dev(sigma)).trace()

def sigma_mag(sigma):
  #Returns the magnitude of a second-rank tensor
  #return np.linalg.norm(sigma)
  return np.sqrt(DblDot(sigma,sigma))

def DblDot(x,y):#Returns the double inner product of two second-rank tensors
  val=0
  for i in range(0,3):
    for j in range(0,3):
      val=val+(x[i][j]*y[i][j])
  return val

def sigma_tau(sigma):
  #return sign(np.sqrt(sigma_J2(sigma)),sigma_J3(sigma))
  return sign(np.sqrt(sigma_J2(sigma)),sigma_J3(sigma))

def get_ps_and_qs(sigmas):
  ps = []
  qs = []
  for sigma in sigmas:
    qs.append(sign(sqrtThree*np.sqrt(sigma_J2(sigma)),sigma_J3(sigma)))
    ps.append(sigma_I1(sigma)/3.0)
  return ps,qs
  
#---------------------------------------------------------------------------------
# Read the internal state variables from the UDA file for a set of times
# using linear interpolation
#---------------------------------------------------------------------------------
def getInternalVariables(uda_path, analytical_times, matID = 0):

  # Get the internal variables
  ev_e, ev_p, capX, times = getAllInternalVariables(uda_path, matID)

  # Create a clean list containing the data as functions of time
  ev_e_list = []
  ev_p_list = []
  capX_list = []
  time_list = []
  #print(analytical_times)
  #print(times)
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times)-1])
  #print("Analytical times = ", an_times_add)
  time_list.append(times[0])
  ev_e_list.append(ev_e[0])
  ev_p_list.append(ev_p[0])
  capX_list.append(capX[0])
  for ii in range(1, len(times)):
    tt_0 = times[ii-1]
    tt_1 = times[ii]
    for jj, ta in enumerate(an_times_add):
      ss = (ta - tt_0)/(tt_1 - tt_0)
      #print("ii = " , ii, " tt_0 = " , tt_0, " tt_1 = ", tt_1, " jj = " , jj, " ta = " , ta, " ss = ", ss )
      if (ss > 0.0 and ss < 1.0):
        #print("ii = " , ii, " tt_0 = " , tt_0, " tt_1 = ", tt_1, " jj = " , jj, " ta = " , ta )
        time_list.append(ta)
        ev_e_list.append((1-ss)*ev_e[ii-1]+ss*ev_e[ii]) 
        ev_p_list.append((1-ss)*ev_p[ii-1]+ss*ev_p[ii]) 
        capX_list.append((1-ss)*capX[ii-1]+ss*capX[ii]) 

  time_list.append(times[len(times)-1])
  ev_e_list.append(ev_e[len(times)-1])
  ev_p_list.append(ev_p[len(times)-1])
  capX_list.append(capX[len(times)-1])
   
  #print("time = ", time_list)
  #print("ev_e = ", ev_e_list)
  #print("ev_p = ", ev_p_list)
  #print("capX = ", capX_list)

  return ev_e_list, ev_p_list, capX_list, time_list, ev_e, ev_p, capX

#---------------------------------------------------------------------------------
# Read the internal state variables from the UDA file for a set of times
# using closest point interpolation
#---------------------------------------------------------------------------------
def find_nearest(times, time):
  idx = np.searchsorted(times, time, side='left')
  if idx > 0 and \
      (idx == len(times) or \
      math.fabs(time - times[idx-1]) < math.fabs(time - times[idx])):
    return idx-1, times[idx-1]
  else:
    return idx, times[idx]

def getInternalVariableSnapshots(uda_path, analytical_times, matID = 0):

  # Get the internal variables
  ev_e, ev_p, capX, times = getAllInternalVariables(uda_path, matID)

  # Create a clean list containing the data as functions of time
  idx_list = []
  time_list = []
  ev_e_list = []
  ev_p_list = []
  capX_list = []
  for ta in analytical_times:
    idx, time = find_nearest(times, ta)
    idx_list.append(idx)
    time_list.append(time)
    ev_e_list.append(ev_e[idx]) 
    ev_p_list.append(ev_p[idx]) 
    capX_list.append(capX[idx]) 

  return idx_list, time_list, ev_e_list, ev_p_list, capX_list, ev_e, ev_p, capX

#---------------------------------------------------------------------------------
# Read the internal state variables from the UDA file:
#---------------------------------------------------------------------------------
def getAllInternalVariables(uda_path, matID):

  # Get the internal variables
  times, ev_e_list = get_pElasticStrainVol(uda_path, matID)
  times_list, ev_p_list = get_pPlasticStrainVol(uda_path, matID)
  capX_times_list, capX_list = get_pCapX(uda_path, matID)

  #print(times_list, ev_e_list, ev_p_list, capX_list)

  return ev_e_list, ev_p_list, capX_list, times_list

def get_pStress(uda_path, matID = 0):
  NAN_FAIL = False
  #Extract stress history
  print("Extracting stress history...")
  args = [partextract_exe, "-mat", str(matID), "-partvar","p.stress",uda_path]
  print(args)
  F_stress = tempfile.TemporaryFile()
  #F_stress = open("./tempStressFileOut.txt","w+")
  #open(os.path.split(uda_path)[0]+'/stressHistory.dat',"w+")
  tmp = sub_proc.Popen(args,stdout=F_stress,stderr=sub_proc.PIPE)
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
    sigmas.append(np.array([[S11,S12,S13],[S21,S22,S23],[S31,S32,S33]]))
    for i in range(3):
      for j in range(3):
        if np.isnan(sigmas[-1][i][j]):
          NAN_FAIL = True
  F_stress.close()
  if NAN_FAIL:
    print("\nERROR: 'nan's found reading in stress. Will not plot correctly")
  return times,sigmas

def get_pDeformationMeasure(uda_path, matID = 0):
  NAN_FAIL = False
  #Extract stress history
  print("Extracting deformation history...")
  args = [partextract_exe, "-mat", str(matID), "-partvar","p.deformationGradient",uda_path]
  F_defMes = tempfile.TemporaryFile()
  #open(os.path.split(uda_path)[0]+'/stressHistory.dat',"w+")
  tmp = sub_proc.Popen(args,stdout=F_defMes,stderr=sub_proc.PIPE)
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
    Fs.append(np.array([[F11,F12,F13],[F21,F22,F23],[F31,F32,F33]]))
    for i in range(3):
      for j in range(3):
        if np.isnan(Fs[-1][i][j]):
          NAN_FAIL = True
  F_defMes.close()
  if NAN_FAIL:
    print("\nERROR: 'nan's found reading in stress. Will not plot correctly")
  return times,Fs

def get_epsilons(uda_path):
  #Assumes no shear strains
  times,Fs = get_pDeformationMeasure(uda_path)
  epsils = []
  for F in Fs:
    epsils.append(np.array([[np.log(F[0][0]),0,0],[0,np.log(F[1][1]),0],[0,0,np.log(F[2][2])]]))
  return times,epsils

#---------------------------------------------------------------------
# Extract capX history
#---------------------------------------------------------------------
def get_pCapX(uda_path, matID = 0):
  print("Extracting capX history...")
  FAIL_NAN = False
  args = [partextract_exe, "-mat", str(matID), "-partvar","p.capX",uda_path]
  print(args)
  F_capX = tempfile.TemporaryFile()
  #open(os.path.split(uda_path)[0]+'/plasticStrainVolHistory.dat',"w+")
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
    capX.append(np.float64(line[4]))
    if np.isnan(capX[-1]):
      FAIL_NAN = True
  F_capX.close()
  if FAIL_NAN:
    print("\ERROR: 'nan' encountered while retrieving p.capX, will not plot correctly."  )
  return times, capX  

#---------------------------------------------------------------------
# Extract plasticStrainVol history
#---------------------------------------------------------------------
def get_pPlasticStrainVol(uda_path, matID = 0):
  FAIL_NAN = False
  #Extract stress history
  print("Extracting plasticStrainVol history...")
  args = [partextract_exe, "-mat", str(matID), "-partvar","p.plasticVolStrain",uda_path]
  print(args)
  F_plasticStrainVol = tempfile.TemporaryFile()
  #open(os.path.split(uda_path)[0]+'/plasticStrainVolHistory.dat',"w+")
  tmp = sub_proc.Popen(args,stdout=F_plasticStrainVol,stderr=sub_proc.PIPE)
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
    print("\ERROR: 'nan' encountered while retrieving p.evp, will not plot correctly."  )
  return times,plasticStrainVol  

#---------------------------------------------------------------------
# Extract elasticStrainVol history
#---------------------------------------------------------------------
def get_pElasticStrainVol(uda_path, matID = 0):
  FAIL_NAN = False
  #Extract elastic strain history
  print("Extracting elasticStrainVol history...")
  args = [partextract_exe, "-mat", str(matID), "-partvar","p.elasticVolStrain",uda_path]
  F_elasticStrainVol = tempfile.TemporaryFile()
  #open(os.path.split(uda_path)[0]+'/elasticStrainVolHistory.dat',"w+")
  tmp = sub_proc.Popen(args,stdout=F_elasticStrainVol,stderr=sub_proc.PIPE)
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
    print("\ERROR: 'nan' encountered while retrieving p.eve, will not plot correctly.")
  return times,elasticStrainVol  

def get_totalStrainVol(uda_path):
  times,plasticStrainVol = get_pPlasticStrainVol(uda_path)
  times,elasticStrainVol = get_pElasticStrainVol(uda_path)
  
  print('num plastic : ',len(plasticStrainVol))
  print('num elastic : ',len(elasticStrainVol))

  totalStrainVol = np.array(plasticStrainVol)+np.array(elasticStrainVol)
  return times,totalStrainVol

def get_defTable(uda_path,working_dir):
  #Determine the defTable file
  try:
    ups_file = os.path.abspath(uda_path)+'/input.xml.orig'
    F = open(ups_file,"r")
  except:
    ups_file = os.path.abspath(uda_path)+'/input.xml'
    F = open(ups_file,"r")
  for line in F:
    if '<prescribed_deformation_file>' in line and '</prescribed_deformation_file>' in line:
      def_file = line.split('<prescribed_deformation_file>')[1].split('</prescribed_deformation_file>')[0].strip()
  F.close()
  #Assumes the input deck and uda share the same parent folder. 
  def_file = working_dir+'/'+def_file
  F = open(def_file,'r')
  times = []
  Fs = []
  for line in F:
    line = line.strip().split()
    times.append(float(line[0]))
    Fs.append(np.array([[float(line[1]),float(line[2]),float(line[3])],
      [float(line[4]),float(line[5]),float(line[6])],
      [float(line[7]),float(line[8]),float(line[9])]]))
  F.close()
  return times,Fs 

def exp_fmt(x,loc):
  tmp = format(x,'1.2e').split('e')
  lead = tmp[0]
  exp = str(int(tmp[1]))
  if exp=='0' and lead=='0.00':
    return r'$\mathbf{0.00}$'
  else:
    if int(exp)<10 and int(exp)>0:
      exp = '+0'+exp
    elif int(exp)>-10 and int(exp)<0:
      exp = '-0'+exp.split('-')[1]
    elif int(exp)>10:
      exp = '+'+exp
    return r'$\mathbf{'+lead+r'\cdot{}10^{'+exp+'}}$' 

#------------------------------------------------------------------------
# Plot x vs y
#------------------------------------------------------------------------
def plotEqShearMeanStress(xs, ys, idx_list, color_list, ev_p_list,
                          compression = 'negative', 
                          Xlims = False, Ylims = False, GRID = True):

  ax1 = plt.subplot(111)
  if (compression == 'positive'):
    xs = list(map(lambda x : -x, xs))
    ys = list(map(lambda y : -y, ys))

  x_data = np.array(xs)
  y_data = np.array(ys)

  max_idx = max(idx_list)
  for jj, color in enumerate(color_list):
    ev_p_str = str(ev_p_list[jj])
    label_str = 'Plastic vol. strain = ' + ev_p_str
    start_idx = idx_list[jj]
    if (start_idx < max_idx):
      end_idx = idx_list[jj+1]
      x_arrow = x_data[start_idx:end_idx-1]
      y_arrow = y_data[start_idx:end_idx-1]
      u_arrow = x_data[start_idx+1:end_idx] - x_arrow
      v_arrow = y_data[start_idx+1:end_idx] - y_arrow
      plt.quiver(x_arrow, y_arrow,
                 u_arrow, v_arrow,
                 scale_units='xy', angles='xy', scale=1.5, width=0.001,
                 headwidth=7, headlength=10,
                 linewidth=0.5, edgecolor='b', color='k')
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
    ax1.set_xlim(Xlims[0],Xlims[1])
  if Ylims:
    ax1.set_ylim(Ylims[0],Ylims[1])
  if GRID:    
    plt.grid(True)  

  return ax1

def eqShear_vs_meanStress(xs, ys, compression = 'negative', 
                          Xlims = False, Ylims = False, 
                          LINE_LABEL = 'Vaango', GRID = True):

  ax1 = plt.subplot(111)
  if (compression == 'negative'):
    x_data = np.array(xs)
    y_data = np.array(ys)
    plt.quiver(x_data[:-1], y_data[:-1],
               x_data[1:]-x_data[:-1], y_data[1:]-y_data[:-1],
               scale_units='xy', angles='xy', scale=1.0, width=0.001,
               headwidth=7, headlength=10,
               linewidth=0.5, edgecolor='b', color='k')
    plt.plot(x_data, y_data,'-r',label=LINE_LABEL)
  else:
    xs = list(map(lambda x : -x, xs))
    ys = list(map(lambda y : -y, ys))
    x_data = np.array(xs)
    y_data = np.array(ys)
    plt.quiver(x_data[:-1], y_data[:-1],
               x_data[1:]-x_data[:-1], y_data[1:]-y_data[:-1],
               scale_units='xy', angles='xy', scale=1.0, width=0.001,
               headwidth=7, headlength=10,
               linewidth=0.5, edgecolor='b', color='k')
    plt.plot(x_data, y_data,'-r',label=LINE_LABEL)
  
  plt.xlabel(str_to_mathbf('Mean Stress, p (Pa)'))
  plt.ylabel(str_to_mathbf('Equivalent Shear Stress, q, (Pa)'))
  
  formatter_int = ticker.FormatStrFormatter('$\mathbf{%g}$')
  formatter_exp = ticker.FuncFormatter(exp_fmt)
  
  #ax1.xaxis.set_major_formatter(formatter_exp)
  #ax1.yaxis.set_major_formatter(formatter_exp)    
  ax1.xaxis.set_major_formatter(formatter_int)
  ax1.yaxis.set_major_formatter(formatter_int)    

  if Xlims:
    ax1.set_xlim(Xlims[0],Xlims[1])
  if Ylims:
    ax1.set_ylim(Ylims[0],Ylims[1])
  if GRID:    
    plt.grid(True)  

  return ax1


def get_rs(nPoints,FSLOPE,PEAKI1,P0,CR):

  kappa = get_kappa(PEAKI1,P0,CR)  
  I1s = np.linspace(PEAKI1,P0,nPoints)
  rs = []
  for I1 in I1s:
    inner_root = (1.0-(pow(kappa-I1,2.0)/pow(kappa-P0,2.0)))
    r = FSLOPE*(I1-PEAKI1)*np.sqrt(2.0*inner_root)
    rs.append(r)
  return I1s,rs
  
def I1_to_zbar(I1s):
  sqrt_3 = np.sqrt(3.0)
  if type(I1s) in [list,np.ndarray]:
    zbars = []
    for I1 in I1s:
      zbars.append(-I1/sqrt_3)
    return zbars
  elif type(I1s) in [int,float,np.float64]:
    return -I1s/sqrt_3
  else:
    print('\nERROR: cannot compute zbar from I1. Invalid type.\n\ttype(I1)\t:\t',type(I1s))
    return None


def defTable_to_J2Solution(def_times, Fs, bulk_mod, shear_mod, tau_yield, num_substeps=1000):
  #Assumes: 
  print('Solving for analytical solution...')
  analytical_epsils = [np.array([[0,0,0],[0,0,0],[0,0,0]])]
  analytical_sigmas = [np.array([[0,0,0],[0,0,0],[0,0,0]])]
  analytical_times = [def_times[0]]
  
  epsils = []
  for F in Fs:
    epsils.append(np.array([[np.log(sum(F[0])),0,0],[0,np.log(sum(F[1])),0],[0,0,np.log(sum(F[2]))]]))
    
  for leg in range(len(def_times)-1):
    t_start = def_times[leg]
    leg_delT = def_times[leg+1]-t_start
    leg_sub_delT = float(leg_delT)/float(num_substeps)
    leg_del_epsil = (epsils[leg+1]-epsils[leg])
    leg_epsil_dot = leg_del_epsil/leg_delT
    for i in range(num_substeps):
      t_now = t_start+float(i)*leg_sub_delT
      analytical_times.append(t_now)
      analytical_sigmas.append(J2VM(leg_epsil_dot,leg_sub_delT,analytical_sigmas[-1],bulk_mod,shear_mod,tau_yield))
      analytical_epsils.append(analytical_epsils[-1]+(leg_epsil_dot*leg_sub_delT))
  analytical_epsils.append(analytical_epsils[-1]+(leg_epsil_dot*leg_sub_delT))   
  analytical_sigmas.append(J2VM(leg_epsil_dot,leg_sub_delT,analytical_sigmas[-1],bulk_mod,shear_mod,tau_yield))
  analytical_times.append(def_times[-1])
  print('Done.')
  return analytical_times,analytical_sigmas,analytical_epsils

def J2_at_Yield(uda_path):
  material_dict = get_yield_surface_data(uda_path)
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
  
  kappa_initial = get_kappa(PEAKI1,P0,CR)  
  I1 = 0
  I1_plus3Pf0 = I1+3.0*Pf0
  if I1_plus3Pf0 >= kappa_initial and I1_plus3Pf0<= PEAKI1:
    J2 = (FSLOPE**2)*((I1-PEAKI1+3.0*Pf0)**2)
  elif I1_plus3Pf0 >= P0 and I1_plus3Pf0 < kappa_initial:
    J2 = ((FSLOPE**2)*((I1-PEAKI1+3.0*Pf0)**2))*(1.0-((I1+CR*FSLOPE*I1-P0-CR*FSLOPE*PEAKI1+3.0*Pf0+3.0*CR*FSLOPE*Pf0)**2/((CR**2)*(FSLOPE**2)*(P0-PEAKI1)**2)))
  else:
    J2 = 0.0  
  return J2
  
def plot_yield_surface(uda_path,PLOT_TYPE='J2_vs_I1'):
  num_points = 500
  material_dict = get_yield_surface_data(uda_path)
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
  kappa_initial = get_kappa(PEAKI1,P0,CR)
  I1s = np.linspace(P0-3.0*Pf0,PEAKI1-3.0*Pf0,num_points)
  #print('Region 1:: ','I1 >= kappa initial-3.0*Pf0 : ',kappa_initial-3.0*Pf0,' ','I1 <= PEAKI1-3*Pf0 : ',PEAKI1-3.0*Pf0)
  #print('Region 2:: ','I1 >= P0-3*Pf0 : ',P0-3.0*Pf0,' ','I1 < kappa_initial-3*Pf0 : ',kappa_initial-3.0*Pf0)
  #print('Region 3:: Not Region 1 or 2')

  #J2 versus I1
  J2s = []
  PLOT = True
  for I1 in I1s:
    I1_plus3Pf0 = I1+3.0*Pf0
    if I1_plus3Pf0 >= kappa_initial and I1_plus3Pf0<= PEAKI1:
      J2 = (FSLOPE**2)*((I1-PEAKI1+3.0*Pf0)**2)
    elif I1_plus3Pf0 >= P0 and I1_plus3Pf0 < kappa_initial:
      Ff = FSLOPE*(PEAKI1-I1_plus3Pf0)  
      fc = np.sqrt(1.0 - ( (kappa_initial-I1_plus3Pf0)/(kappa_initial-P0) )**2)
      J2 = (Ff*fc)**2
    else:
      J2 = 0.0
    J2s.append(J2)
  
  if PLOT_TYPE == 'J2_vs_I1':
    xs = np.array(I1s)
    ys = np.array(J2s)*1.0e-12  
  elif PLOT_TYPE == 'sqrtJ2_vs_I1':
    xs = np.array(I1s)
    ys = np.sqrt(np.array(J2s))
  elif PLOT_TYPE == 'r_vs_z':
    xs = np.array(I1s)/np.sqrt(3.0)
    ys = np.sqrt(2.0*np.array(J2s))
  elif PLOT_TYPE == 'q_vs_I1':
    xs = np.array(I1s)
    ys = np.sqrt(3.0*np.array(J2s))
  elif PLOT_TYPE == 'q_vs_p':
    xs = np.array(I1s)/3.0
    ys = np.sqrt(3.0*np.array(J2s))
  else:
    PLOT = False
    print('\nError: invalid plot type specified for initial yield surface plot.\n\tPLOT_TYPE:',PLOT_TYPE)
  if PLOT:
    plt.plot(xs,ys,'--k',linewidth=lineWidth+1,label='Initial Yield Surface')
    plt.plot(xs,-ys,'--k',linewidth=lineWidth+1)  

def test_yield_surface(uda_path):
  plot_yield_surface_2(uda_path,'J2_vs_I1')
  plt.show()
  plot_yield_surface_2(uda_path,'sqrtJ2_vs_I1')
  plt.show()
  plot_yield_surface_2(uda_path,'r_vs_z')
  plt.show()
  plot_yield_surface_2(uda_path,'q_vs_I1')
  plt.show()
  plot_yield_surface_2(uda_path,'q_vs_p')
  plt.show()
  
def plot_yield_surface_updated(uda_path, ev_p, zeta=0.0, PLOT_TYPE='J2_vs_I1', COLOR='#0000ff'):
  num_points = 500
  material_dict = get_yield_surface_data(uda_path)
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
  a2 = (FSLOPE-YSLOPE)/(STREN-YSLOPE*PEAKI1) 
  a3 = (STREN-YSLOPE*PEAKI1)*math.exp(-a2*PEAKI1)
  a4 = YSLOPE
  phi_i = 1.0 - math.exp(-P3);      
  Km = B0 + B1;                       
  Kf = fluid_B0;                           
  C1 = Kf*(1.0 - phi_i) + Km*(phi_i);   
  print('Kf = ', Kf, ' Km = ', Km)
  ev0 = 0.0;
  if (Kf != 0):
    ev0 = C1*Pf0/(Kf*Km)

  # compute capX
  cap_X = computeCapX(ev_p, P3, P0, P1, Kf, Km, ev0, C1, phi_i, B0, B1, B2, B3, B4)
  print("capX = ", cap_X)

  # Compute kappa
  kappa = PEAKI1 - CR*(PEAKI1 - cap_X)
  kappa_initial = PEAKI1 - CR*(PEAKI1 - P0)
  print("kappa = ", kappa, " kappa_initial = ", kappa_initial)

  #I1s = np.linspace(P0-3.0*Pf0,PEAKI1-3.0*Pf0,num_points)
  I1s = np.linspace(0.999*cap_X - 3.0*Pf0, PEAKI1-3.0*Pf0, num_points)
  #print('Region 1:: ','I1 >= kappa initial-3.0*Pf0 : ',kappa_initial-3.0*Pf0,' ','I1 <= PEAKI1-3*Pf0 : ',PEAKI1-3.0*Pf0)
  #print('Region 2:: ','I1 >= P0-3*Pf0 : ',P0-3.0*Pf0,' ','I1 < kappa_initial-3*Pf0 : ',kappa_initial-3.0*Pf0)
  #print('Region 3:: Not Region 1 or 2')

  #J2 versus I1
  J2s = []
  PLOT = True
  for I1 in I1s:
    I1f = I1+3.0*Pf0

    # Compute F_f
    I1_minus_zeta = I1f - zeta;
    Ff = a1 - a3*math.exp(a2*I1_minus_zeta) - a4*(I1_minus_zeta)
    Ff_sq = Ff**2

    # Compute Fc
    Fc_sq = 1.0
    if (I1_minus_zeta < kappa) and (cap_X < I1_minus_zeta):
      ratio = (kappa - I1_minus_zeta)/(kappa - cap_X)
      Fc_sq = 1.0 - ratio**2

    # Compute J2
    J2 = Ff_sq*Fc_sq
    J2s.append(J2)
  
  print(len(I1s), len(J2s))
  if PLOT_TYPE == 'J2_vs_I1':
    xs = np.array(I1s)
    ys = np.array(J2s)*1.0e-12  
  elif PLOT_TYPE == 'sqrtJ2_vs_I1':
    xs = np.array(I1s)
    ys = np.sqrt(np.array(J2s))
  elif PLOT_TYPE == 'r_vs_z':
    xs = np.array(I1s)/np.sqrt(3.0)
    ys = np.sqrt(2.0*np.array(J2s))
  elif PLOT_TYPE == 'q_vs_I1':
    xs = np.array(I1s)
    ys = np.sqrt(3.0*np.array(J2s))
  elif PLOT_TYPE == 'q_vs_p':
    xs = np.array(I1s)/3.0
    ys = np.sqrt(3.0*np.array(J2s))
  else:
    PLOT = False
    print('\nError: invalid plot type specified for updated yield surface plot.\n\tPLOT_TYPE:',PLOT_TYPE)
  if PLOT:
    print(xs.shape, ys.shape)
    ev_p_str = str(ev_p)
    label_str = 'Yield surf. evp = ' + ev_p_str
    line1 = plt.plot(xs,ys,'--b',linewidth=lineWidth+1,label=label_str)
    line2 = plt.plot(xs,-ys,'--b',linewidth=lineWidth+1)  
    plt.setp(line1, color=COLOR)
    plt.setp(line2, color=COLOR)

def test_yield_surface(uda_path):
  plot_yield_surface_2(uda_path,'J2_vs_I1')
  plt.show()
  plot_yield_surface_2(uda_path,'sqrtJ2_vs_I1')
  plt.show()
  plot_yield_surface_2(uda_path,'r_vs_z')
  plt.show()
  plot_yield_surface_2(uda_path,'q_vs_I1')
  plt.show()
  plot_yield_surface_2(uda_path,'q_vs_p')
  plt.show()
  
###
# Compute elastic
###   
def computeElastic(I1, ev_p, Kf, Km, ev0, C1, phi_i, B0, B1, B2, B3, B4, G0):

  bulk = B0
  shear = G0

  if (ev_p <= 0.0):

    if (ev_p < 0.0):
      bulk = bulk - B3*math.exp(B4/ev_p)

    if (I1 < 0.0):
      bulk = bulk + B1*math.exp(B2/I1)

  if (ev_p <= ev0 and Kf != 0.0):
    Kd = B0;
    if (ev_p < 0.0):
      Kd = B0 - B3*math.exp(B4/ev_p)
      C2 = math.exp(ev_p*Km/C1)*phi_i
      phi = C2/(-math.exp(ev_p*Kf/C1)*(phi_i-1.0) + C2)
      tmp = 1.0 - Kd/Km
      bulk = Kd + tmp**2/((tmp - phi)/Km + (1.0/Kf-1.0/Km)*phi)

  return bulk, shear

#---------------------------------------------------------------------------------
# Read the elastic model/yield function table from json file
#---------------------------------------------------------------------------------
def getJSONTable(json_file):

  with open(json_file) as f:
    data = json.load(f)
  return data['Vaango_tabular_data']['Data']

#-----------------------------------------------------------------------------------
# Read the simulation stress data and compute p,q (converts to Pa)
#-----------------------------------------------------------------------------------
def readSimStressData(uda_path, matID = 0):

  NAN_FAIL = False

  #Extract stress history
  print("Extracting stress history...")
  args = [partextract_exe, "-mat", str(matID), "-partvar","p.stress",uda_path]
  print(args)
  F_stress = tempfile.TemporaryFile()
  tmp = sub_proc.Popen(args,stdout=F_stress,stderr=sub_proc.PIPE)
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
    #print("Line = ", line)
    first_word = line.decode('utf8').partition(' ')[0]
    if first_word == "Error":
      print("**ERROR** partextract failed to read the stress history data.")
      print("          It generated the following message:")
      print(line)
      sys.exit("Stopping program.")

    line = line.strip().split()
    S11 = np.float64(line[4])
    S12 = np.float64(line[5])
    S13 = np.float64(line[6])
    S21 = np.float64(line[7])
    S22 = np.float64(line[8])
    S23 = np.float64(line[9])
    S31 = np.float64(line[10])
    S32 = np.float64(line[11])
    S33 = np.float64(line[12])

    sigma_a = S11
    sigma_r = S22
    sigma_ar = S12
    
    time_sim.append(float(line[0]))
    sigma_sim.append(np.array([[S11,S12,S13],[S21,S22,S23],[S31,S32,S33]]))
    #print("stress = ", sigma_sim)
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
  J3_sim = []
  for sigma in sigma_sim:
    I1_sim.append(sigma_I1(sigma))
    J2_sim.append(sigma_J2(sigma))
    J3_sim.append(sigma_J3(sigma))

  # Compute p, q
  pp_sim = list(map(lambda I1 : I1/3, I1_sim))
  qq_sim = list(map(lambda J2, J3 : np.sign(J3)*np.sqrt(3*J2), J2_sim, J3_sim))

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
def plotSimDataSigmaEps(fig, time_snapshots, time_sim, pp_sim, ev_e_sim, ev_p_sim,
                        compression = 'negative'):

  # Get snapshots from data
  time_snap, pp_snap = getDataTimeSnapshots(time_snapshots, time_sim, pp_sim)

  # Activate the figure
  plt.figure(fig.number)

  ev_sim = list(map(lambda ev_e, ev_p : ev_e + ev_p, ev_e_sim ,ev_p_sim))

  # Plot sigma_a vs. time
  if (compression == 'positive'):
    ev_data = list(map(lambda p : -p, ev_sim))
    pp_data = list(map(lambda p : -p, pp_sim))
    plt.plot(ev_data, pp_data, '--r', label='Simulation')
  else:
    plt.plot(ev_sim, pp_sim, '--r', label='Simulation')

  return time_snap, pp_snap

#-----------------------------------------------------------------------------------
# Plot the sim data as a function of time (assume Pa)
#-----------------------------------------------------------------------------------
def plotSimDataSigmaTime(fig, time_snapshots, time_sim, sigma_a_sim, sigma_r_sim, sigma_ar_sim,
                         labelxx = '$\sigma_a$ (sim)', labelyy = '$\sigma_r$ (sim)',
                         labelxy = '$\sigma_{ar}$ (sim)',
                         compression = 'negative'):

  # Get snapshots from data
  time_snap, sigma_a_snap = getDataTimeSnapshots(time_snapshots, time_sim, sigma_a_sim)
  time_snap, sigma_r_snap = getDataTimeSnapshots(time_snapshots, time_sim, sigma_r_sim)
  time_snap, sigma_ar_snap = getDataTimeSnapshots(time_snapshots, time_sim, sigma_ar_sim)

  # Activate the figure
  plt.figure(fig.number)

  # Plot sigma_a vs. time
  if (compression == 'positive'):
    sigma_a_data = list(map(lambda p : -p, sigma_a_sim))
    sigma_r_data = list(map(lambda p : -p, sigma_r_sim))
    sigma_ar_data = list(map(lambda p : -p, sigma_ar_sim))
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
    plt_color = cm.Paired(float(ii)/len(time_snap))
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
def plotSimDataPQTime(fig, time_snapshots, time_sim, p_sim, q_sim, compression='negative'):

  # Get snapshots from sim data
  time_snap, p_snap = getDataTimeSnapshots(time_snapshots, time_sim, p_sim)
  time_snap, q_snap = getDataTimeSnapshots(time_snapshots, time_sim, q_sim)

  # Activate the figure
  plt.figure(fig.number)

  # Plot sigma_a vs. time
  if (compression == 'positive'):
    p_data = list(map(lambda p : -p, p_sim))
    q_data = list(map(lambda q : -q, q_sim))
    plt.plot(time_sim, p_data, '--r', label='$p$ (sim)')
    plt.plot(time_sim, q_data, '--b', label='$q$ (sim)')
  else:
    plt.plot(time_sim, p_sim, '--r', label='$p$ (sim)')
    plt.plot(time_sim, q_sim, '--b', label='$q$ (sim)')

  # Plot filled circles at time snapshots
  for ii in range(0, len(time_snap)):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii)/len(time_snap))
    if (compression == 'positive'):
      plt.plot(time_snap[ii], -p_snap[ii], 'o', color=plt_color) 
      plt.plot(time_snap[ii], -q_snap[ii], 'o', color=plt_color) 
    else:
      plt.plot(time_snap[ii], p_snap[ii], 'o', color=plt_color) 
      plt.plot(time_snap[ii], q_snap[ii], 'o', color=plt_color) 

  return time_snap, p_snap, q_snap

