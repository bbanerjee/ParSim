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

from TabularYieldSurfaceUtils import *

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
  plt.savefig(name+'.png', bbox_inches=0, dpi=plt.rcParams['figure.dpi']) 

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
def getInternalVariables(uda_path, analytical_times):

  # Get the internal variables
  ev_e, ev_p, times = getAllInternalVariables(uda_path)

  # Create a clean list containing the data as functions of time
  ev_e_list = []
  ev_p_list = []
  time_list = []
  #print(analytical_times)
  #print(times)
  an_times_add = list(analytical_times)
  an_times_add.append(times[len(times)-1])
  #print("Analytical times = ", an_times_add)
  time_list.append(times[0])
  ev_e_list.append(ev_e[0])
  ev_p_list.append(ev_p[0])
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
    #    break
    #else:
    #  continue
  time_list.append(times[len(times)-1])
  ev_e_list.append(ev_e[len(times)-1])
  ev_p_list.append(ev_p[len(times)-1])
   
  #print("time = ", time_list)
  #print("ev_e = ", ev_e_list)
  #print("ev_p = ", ev_p_list)

  return ev_e_list, ev_p_list, time_list

#---------------------------------------------------------------------------------
# Read the internal state variables from the UDA file:
#---------------------------------------------------------------------------------
def getAllInternalVariables(uda_path):

  # Get the internal variables
  times, ev_e_list = get_pElasticStrainVol(uda_path)
  times, ev_p_list = get_pPlasticStrainVol(uda_path)

  return ev_e_list, ev_p_list, times

def get_pStress(uda_path):
  NAN_FAIL = False
  #Extract stress history
  print("Extracting stress history...")
  args = [partextract_exe, "-partvar","p.stress",uda_path]
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

def get_pDeformationMeasure(uda_path):
  NAN_FAIL = False
  #Extract stress history
  print("Extracting deformation history...")
  args = [partextract_exe, "-partvar","p.deformationGradient",uda_path]
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

def get_pPlasticStrainVol(uda_path):
  FAIL_NAN = False
  #Extract stress history
  print("Extracting plasticStrainVol history...")
  args = [partextract_exe, "-partvar","p.plasticVolStrain",uda_path]
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

def get_pElasticStrainVol(uda_path):
  FAIL_NAN = False
  #Extract elastic strain history
  print("Extracting elasticStrainVol history...")
  args = [partextract_exe, "-partvar","p.elasticVolStrain",uda_path]
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

def eqShear_vs_meanStress(xs,ys,Xlims=False,Ylims=False,LINE_LABEL='Vaango',GRID=True):

  ax1 = plt.subplot(111)
  plt.plot(np.array(xs),np.array(ys),'-r',label=LINE_LABEL)
  
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

### ----------
#  Test Methods Below
### ----------
def test02_postProc(uda_path,save_path,**kwargs):
  #Extract stress history
  print("Post Processing Test: 02 - Vertex Treatment")
  times,sigmas = get_pStress(uda_path)
  ps_unscaled,qs_unscaled = get_ps_and_qs(sigmas)
  ps = []
  qs = []
  for val in ps_unscaled:
    ps.append(val)

  for val in qs_unscaled:
    qs.append(val)
 
  Sxx = []
  Syy = []
  Szz = []
  for sigma in sigmas:
    Sxx.append(sigma[0][0])
    Syy.append(sigma[1][1])
    Szz.append(sigma[2][2])

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
  Sxx_tick_int = (Sxx_max - Sxx_min)/4;
  Syy_tick_int = (Syy_max - Syy_min)/4;
  Szz_tick_int = (Szz_max - Szz_min)/4;
  Sxx_ticks = [round(Sxx_max+Sxx_tick_int,1), round(Sxx_max,1), 
               round(Sxx_max-Sxx_tick_int,1), round(Sxx_max-2*Sxx_tick_int,1), 
               round(Sxx_max-3*Sxx_tick_int,1), round(Sxx_max-4*Sxx_tick_int,1)] 
  Syy_ticks = [round(Syy_max+Syy_tick_int,1), round(Syy_max,1), 
               round(Syy_max-Syy_tick_int,1), round(Syy_max-2*Syy_tick_int,1), 
               round(Syy_max-3*Syy_tick_int,1), round(Syy_max-4*Syy_tick_int,1)] 
  Szz_ticks = [round(Szz_max+Szz_tick_int,1), round(Szz_max,1), 
               round(Szz_max-Szz_tick_int,1), round(Szz_max-2*Szz_tick_int,1), 
               round(Szz_max-3*Szz_tick_int,1), round(Szz_max-4*Szz_tick_int,1)] 
  
  #Analytical Solutions
  #Drucker-Prager constants
  r0 = 50.0
  z0 = 50.0*sqrtThree
  #Solution From Brannon Leelavanichkul paper
  analytical_times = [0.0,1.0,threeHalves,2.0,5.0/2.0,3.0]
  analytical_S11 = np.array([0,-850.0/3.0,(-50.0/3.0)*(9.0+4.0*np.sqrt(6.0)),(-50.0/3.0)*(9.0+4.0*np.sqrt(6.0)),(50.0/3.0)*(2.0*np.sqrt(6)-3.0),160.0*np.sqrt(twoThirds)-110.0])
  analytical_S22 = np.array([0,-850.0/3.0,(50.0/3.0)*(2.0*np.sqrt(6.0)-9.0),(50.0/3.0)*(2.0*np.sqrt(6.0)-9.0),(-50.0/3.0)*(3.0+np.sqrt(6.0)),(-10.0/3.0)*(33.0+8.0*np.sqrt(6.0))])
  analytical_S33 = np.array([0,-850.0/3.0,(50.0/3.0)*(2.0*np.sqrt(6.0)-9.0),(50.0/3.0)*(2.0*np.sqrt(6.0)-9.0),(-50.0/3.0)*(3.0+np.sqrt(6.0)),(-10.0/3.0)*(33.0+8.0*np.sqrt(6.0))])
  analytical_mean = (analytical_S11+analytical_S22+analytical_S33)/3.0
  analytical_I1 = analytical_S11+analytical_S22+analytical_S33
  tmp = (1.0/3.0)*analytical_I1
  analytical_s1 = analytical_S11-tmp
  analytical_s2 = analytical_S22-tmp
  analytical_s3 = analytical_S33-tmp
  analytical_J2 = (1.0/2.0)*(pow(analytical_s1,2)+pow(analytical_s2,2)+pow(analytical_s3,2))
  analytical_J3 = (1.0/3.0)*(pow(analytical_s1,3)+pow(analytical_s2,3)+pow(analytical_s3,3))
  analytical_z = analytical_I1/sqrtThree
  analytical_q = []
  for idx,J2 in enumerate(analytical_J2):
    J3 = analytical_J3[idx]
    analytical_q.append(sign(sqrtThree*np.sqrt(J2),J3))
  #Drucker-Prager yield surface
  yield_zs = np.array([z0,min(analytical_z)])
  yield_rs = r0/z0*((get_yield_surface_data(uda_path)['PEAKI1']/sqrtThree)-yield_zs)
  yield_ps = yield_zs*(sqrtThree/3.0)
  yield_qs = yield_rs*np.sqrt(threeHalves)

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
  an_times_add.append(times[len(times)-1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4  and jj > 0):
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
  material_dict = get_yield_surface_data(uda_path)  
  param_text = material_dict['material string']
  plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')  
  #plt.plot(analytical_mean,analytical_q,'-g',linewidth=lineWidth+1,label='Analytical')
  #plt.plot(yield_ps,yield_qs,'--k',linewidth=lineWidth+2,label='Yield surface')
  #plt.plot(yield_ps,-yield_qs,'--k',linewidth=lineWidth+2)
  #eqShear_vs_meanStress(ps,qs,(-300,60),(-300,300))  
  #eqShear_vs_meanStress(ps,qs,(ps_min, ps_max),(qs_min, -qs_min))  

  eqShear_vs_meanStress(ps,qs)  
  
  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('TabularTest 02:\nVertex Treatment (plot a)')  
  #plt.legend()
  savePNG(save_path+'/Test02_verificationPlot_a','1280x960')
  
  ##Plot b
  plt.figure(2)
  plt.clf()
  plt.subplots_adjust(right=0.75)
  param_text = material_dict['material string']
  plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')   
  endT = max(times)  
  #Sigma zz
  ax3 = plt.subplot(313)
  plt.plot(analytical_times,analytical_S33,'-g',linewidth=lineWidth+2)
  plt.plot(times,np.array(Szz),'-r')  
  #Add Yield Surface
  #Add Analytical
  plt.xlabel(str_to_mathbf('Time (s)'))  
  plt.ylabel(str_to_mathbf('\sigma_{zz} (Pa)')) 
  ax3.yaxis.set_major_formatter(formatter)  
  ax3.set_xlim(0,endT)
  #ax3.set_ylim(-300,100)
  #ax3.set_yticks([-300,-200,-100,0,100])
  plt.grid(True)
  #Sigma xx
  ax1 = plt.subplot(311,sharex=ax3)
  plt.plot(analytical_times,analytical_S11,'-g',linewidth=lineWidth+2,label='Analytical')
  plt.plot(times,np.array(Sxx),'-r',label='Uintah')  
  #Add Yield Surface
  #Add Analytical  
  plt.legend()
  plt.setp(ax1.get_xticklabels(), visible=False)
  plt.ylabel(str_to_mathbf('\sigma_{xx} (Pa)'))  
  plt.title('TabularTest 02:\nVertex Treatment (plot b)')
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)    
  ax1.set_xlim(0,endT)
  #ax1.set_ylim(-400,100)
  #ax1.set_yticks([-400,-300,-200,-100,0,100])
  plt.grid(True)
  #Sigma yy
  ax2 = plt.subplot(312,sharex=ax3)
  plt.plot(analytical_times,analytical_S22,'-g',linewidth=lineWidth+2)
  plt.plot(times,np.array(Syy),'-r')  
  #Add Yield Surface
  #Add Analytical 
  plt.setp(ax2.get_xticklabels(), visible=False)
  plt.ylabel(str_to_mathbf('\sigma_{yy} (Pa)'))  
  ax2.yaxis.set_major_formatter(formatter)    
  ax2.set_xlim(0,endT)
  #ax2.set_ylim(-300,100)
  #ax2.set_yticks([-300,-200,-100,0,100])  
  plt.grid(True)
  savePNG(save_path+'/Test02_verificationPlot_b','1280x960')   
  
  if SHOW_ON_MAKE:
    plt.show()
  
def test03_postProc(uda_path,save_path,**kwargs):    
  #Extract stress history
  print("Post Processing Test: 03 - Uniaxial Strain Without Hardening")
  times, sigmas = get_pStress(uda_path)
  ps_unscaled, qs_unscaled = get_ps_and_qs(sigmas) 
  material_dict = get_yield_surface_data(uda_path)
  PEAKI1 = material_dict['PEAKI1']
  J2Yield = J2_at_Yield(uda_path)
  q_yield = np.sqrt(3.0*J2Yield)
  
  # Scale the data
  ps = []
  qs = []
  for val in ps_unscaled:
    ps.append(val)

  for val in qs_unscaled:
    qs.append(val)

  PEAKI1 = PEAKI1
  q_yield = q_yield
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
  
  analytical_times = [0.0,1.0,2.0]
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
  an_times_add.append(times[len(times)-1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4  and jj > 0):
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
  material_dict = get_yield_surface_data(uda_path)  
  param_text = material_dict['material string']
  plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')   
  #eqShear_vs_meanStress(ps,qs,Xlims,Ylims,)
  eqShear_vs_meanStress(ps,qs)
  plt.title('TabularTest 03:\nUniaxial Strain Without Hardening')

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.plot(Xlims,(q_yield,q_yield),'--k',linewidth=lineWidth+1,label='Initial yield surface')
  plt.plot(Xlims,(-q_yield,-q_yield),'--k',linewidth=lineWidth+1)
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)   
  #plt.legend()
  savePNG(save_path+'/Test03_verificationPlot','1280x960')
  if SHOW_ON_MAKE:
    plt.show()

def test04_postProc(uda_path,save_path,**kwargs):
  #Extract stress history
  print("Post Processing Test: 04 - Curved Yield Surface")
  times,sigmas = get_pStress(uda_path)
  ps_unscaled,qs_unscaled = get_ps_and_qs(sigmas)

  # Scale the data
  ps = []
  qs = []
  for val in ps_unscaled:
    ps.append(val)

  for val in qs_unscaled:
    qs.append(val)

  # Find min/max values
  ps_min = min(ps)
  ps_max = max(ps)
  qs_min = min(qs)
  qs_max = max(qs)
  print("ps_min = ", ps_min)
  print("ps_max = ", ps_max)
  print("qs_min = ", qs_min)
  print("qs_max = ", qs_max)
  
  analytical_times = [0.0,1.0]
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
  an_times_add.append(times[len(times)-1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4  and jj > 0):
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
  material_dict = get_yield_surface_data(uda_path)  
  param_text = material_dict['material string']
  plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')  
  eqShear_vs_meanStress(ps,qs)

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('TabularTest 04:\nCurved Yield Surface')
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)   
  #Add Analytical
  #plt.legend()
  savePNG(save_path+'/Test04_verificationPlot','1280x960')
  if SHOW_ON_MAKE:
    plt.show()
  
def test05_postProc(uda_path,save_path,**kwargs):
  #Extract stress history
  print("Post Processing Test: 05 - Hydrostatic Compression Fixed Cap")
  times,sigmas = get_pStress(uda_path)
  ps_unscaled,qs_unscaled = get_ps_and_qs(sigmas)

  # Scale the data
  ps = []
  qs = []
  for val in ps_unscaled:
    ps.append(val)

  for val in qs_unscaled:
    qs.append(val)

  # Find min/max values
  ps_min = min(ps)
  ps_max = max(ps)
  qs_min = min(qs)
  qs_max = max(qs)
  print("ps_min = ", ps_min)
  print("ps_max = ", ps_max)
  print("qs_min = ", qs_min)
  print("qs_max = ", qs_max)
  
  analytical_times = [0.0,1.0]
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
  an_times_add.append(times[len(times)-1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 5.0e-4  and jj > 0):
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
  material_dict = get_yield_surface_data(uda_path)  
  param_text = material_dict['material string']
  plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')  
  eqShear_vs_meanStress(ps,qs)

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('TabularTest 05:\nHydrostatic Compression Fixed Cap')
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter) 
  #Add Analytical
  #plt.legend()
  savePNG(save_path+'/Test05_verificationPlot','1280x960')
  if SHOW_ON_MAKE:
    plt.show()
  
def test08_postProc(uda_path,save_path,**kwargs):
  #Extract stress history
  print("Post Processing Test: 08 - Loading/Unloading")
  times,sigmas = get_pStress(uda_path)
  I1s = []
  ps = []
  for sigma in sigmas:
    I1s.append(sigma_I1(sigma))
    ps.append(sigma_I1(sigma)/3.0)

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
  an_times_add.append(times[len(times)-1])
  for ii, tt in enumerate(times):
    for jj, ta in enumerate(an_times_add):
      if (math.fabs(tt - ta) < 1.0e-4  and jj > 0):
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

  totalStrainVol = np.array(elasticStrainVol)+np.array(plasticStrainVol)
  material_dict = get_yield_surface_data(uda_path)
  P3 = material_dict['P3']
  porosity = 1-np.exp(-(P3+np.array(plasticStrainVol)))

  ###PLOTTING
  int_formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  exp_formatter = ticker.FuncFormatter(exp_fmt)
  ##Plot a
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75,left=0.15)  
  param_text = material_dict['material string']
  plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')  
  print(len(times))
  print(len(ps))
  ax1=eqShear_vs_meanStress(times,-np.array(ps))

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('TabularTest 08:\nLoading/Unloading (plot a)')  
  plt.ylabel(str_to_mathbf('Pressure (Pa)'))
  plt.xlabel(str_to_mathbf('Time (s)'))
  ax1.xaxis.set_major_formatter(int_formatter)
  ax1.yaxis.set_major_formatter(exp_formatter)
  ax1.tick_params(axis='both',labelsize='small')
  savePNG(save_path+'/Test08_verificationPlot_a','1280x960')  
  
  ##Plot b
  plt.figure(2)
  plt.clf()
  ax2 = plt.subplot(111)
  plt.subplots_adjust(right=0.75,left=0.15)  
  param_text = material_dict['material string']
  plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')  
  ax1=eqShear_vs_meanStress(times,totalStrainVol)

  # Plot yield surfaces
  plotPQYieldSurfaceSim(uda_path, analytical_times)

  plt.title('TabularTest 08:\nLoading/Unloading (plot b)')  
  plt.ylabel(str_to_mathbf('Total Volumetric Strain, \epsilon_{v}'))
  plt.xlabel(str_to_mathbf('Time (s)'))
  ax2.xaxis.set_major_formatter(int_formatter)
  ax2.yaxis.set_major_formatter(int_formatter)  
  ax2.tick_params(axis='both',labelsize='small')
  savePNG(save_path+'/Test08_verificationPlot_b','1280x960')

  ##Plot c
  I1lims = (I1s_min, I1s_max)  
  plt.figure(3)
  plt.clf()
  ax3 = plt.subplot(111)
  plt.subplots_adjust(right=0.75,left=0.15) 
  param_text = material_dict['material string']
  plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')  
  eqShear_vs_meanStress(I1s,porosity,I1lims,(0,1.25))
  plt.title('TabularTest 08:\nLoading/Unloading (plot c)')
  plt.ylabel(str_to_mathbf('Porosity'))
  plt.xlabel(str_to_mathbf('I_{1}:first invariant of stress tensor (Pa)'))
  plot_crush_curve(uda_path,I1lims)
  #ax1.set_xticks([-9000,-7000,-5000,-3000,-1000,0])
  #ax3.set_xticks([-10000,-7500,-5000,-2500,0,1000])
  ax3.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
  ax3.xaxis.set_major_formatter(exp_formatter)
  ax3.yaxis.set_major_formatter(int_formatter)
  ax3.tick_params(axis='both',labelsize='small')
  plt.legend()
  savePNG(save_path+'/Test08_verificationPlot_c','1280x960')

  if SHOW_ON_MAKE:
    plt.show()   
  
def test10_postProc(uda_path,save_path,**kwargs):
  
  BIG_FIGURE = True
  if 'WORKING_PATH' in kwargs:
    working_dir = kwargs['WORKING_PATH']
  else:
    print('\nERROR: need working directory to post process this problem')
    return

  #Extract stress history
  print("Post Processing Test: 10 - Transient Stress Eigenvalues with Constant Eigenvectors")
  times,sigmas = get_pStress(uda_path)
  Sxx = []
  Syy = []
  Szz = []
  for sigma in sigmas:
    Sxx.append(sigma[0][0])
    Syy.append(sigma[1][1])
    Szz.append(sigma[2][2])  

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
  Sxx_tick_int = (Sxx_max - Sxx_min)/4;
  Syy_tick_int = (Syy_max - Syy_min)/4;
  Szz_tick_int = (Szz_max - Szz_min)/4;
  Sxx_ticks = [round(Sxx_max+Sxx_tick_int,1), round(Sxx_max,1), 
               round(Sxx_max-Sxx_tick_int,1), round(Sxx_max-2*Sxx_tick_int,1), 
               round(Sxx_max-3*Sxx_tick_int,1), round(Sxx_max-4*Sxx_tick_int,1)] 
  Syy_ticks = [round(Syy_max+Syy_tick_int,1), round(Syy_max,1), 
               round(Syy_max-Syy_tick_int,1), round(Syy_max-2*Syy_tick_int,1), 
               round(Syy_max-3*Syy_tick_int,1), round(Syy_max-4*Syy_tick_int,1)] 
  Szz_ticks = [round(Szz_max+Szz_tick_int,1), round(Szz_max,1), 
               round(Szz_max-Szz_tick_int,1), round(Szz_max-2*Szz_tick_int,1), 
               round(Szz_max-3*Szz_tick_int,1), round(Szz_max-4*Szz_tick_int,1)] 
      
  #Analytical solution
  material_dict = get_yield_surface_data(uda_path)
  def_times,Fs = get_defTable(uda_path,working_dir)
  tau_yield = material_dict['PEAKI1']
  bulk_mod = material_dict['B0']
  shear_mod = material_dict['G0']
    
  analytical_times,analytical_sigmas,epsils=defTable_to_J2Solution(def_times,Fs,bulk_mod,shear_mod,tau_yield,num_substeps=10)
   
  analytical_Sxx = []
  analytical_Syy = []
  analytical_Szz = []
  for sigma in analytical_sigmas:
    analytical_Sxx.append(sigma[0][0])
    analytical_Syy.append(sigma[1][1])
    analytical_Szz.append(sigma[2][2])

  ###PLOTTING
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  if BIG_FIGURE:
    plt.subplots_adjust(right=0.75)
    param_text = material_dict['material string']
    plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')
  else:
    plt.subplots_adjust(left=0.15,top=0.96,bottom=0.15,right=0.96)       
   
  #analytical solution
  plt.plot(analytical_times,np.array(analytical_Sxx),':r',linewidth=lineWidth+2,label=str_to_mathbf('Analytical \sigma_{xx}'))  
  plt.plot(analytical_times,np.array(analytical_Syy),'--g',linewidth=lineWidth+2,label=str_to_mathbf('Analytical \sigma_{yy}'))     
  plt.plot(analytical_times,np.array(analytical_Szz),'-.b',linewidth=lineWidth+2,label=str_to_mathbf('Analytical \sigma_{zz}'))    

  #simulation results
  plt.plot(times,np.array(Sxx),'-r',label=str_to_mathbf('Uintah \sigma_{xx}'))       
  plt.plot(times,np.array(Syy),'-g',label=str_to_mathbf('Uintah \sigma_{yy}'))   
  plt.plot(times,np.array(Szz),'-b',label=str_to_mathbf('Uintah \sigma_{zz}'))      
    
  ax1.set_xlim(0,2.25)   
  formatter_int = ticker.FormatStrFormatter('$\mathbf{%g}$')
  ax1.xaxis.set_major_formatter(formatter_int)
  ax1.yaxis.set_major_formatter(formatter_int)     
  #labels
  plt.grid(True)         
  plt.xlabel(str_to_mathbf('Time (s)'))
  plt.ylabel(str_to_mathbf('Stress (Pa)'))
  if BIG_FIGURE:
    plt.legend(loc='upper right', bbox_to_anchor=(1.38,1.12))
    plt.title('TabularTest 10:\nTransient Stress Eigenvalues with Constant Eigenvectors') 
    savePNG(save_path+'/Test10_verificationPlot','1280x960')
  else:
    tmp = plt.rcParams['legend.fontsize']
    plt.rcParams['legend.fontsize']='x-small'
    plt.legend(loc=7)
    savePNG(save_path+'/Test10_verificationPlot','640x480')
    plt.rcParams['legend.fontsize']=tmp

  if SHOW_ON_MAKE:
    plt.show()
  
def test11_postProc(uda_path,save_path,**kwargs):
  if 'WORKING_PATH' in kwargs:
    working_dir = kwargs['WORKING_PATH']
  else:
    print('\nERROR: need working directory to post process this problem'  )
  
  #Extract stress and strain history
  print("Post Processing Test: 11 - Uniaxial Strain J2 Plasticity")
  times,sigmas = get_pStress(uda_path)
  times,epsils = get_epsilons(uda_path)
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
    Sxx.append(sigma[0][0])
    Syy.append(sigma[1][1])
    Szz.append(sigma[2][2])      
      
  #Analytical solution
  material_dict = get_yield_surface_data(uda_path)
  def_times,Fs = get_defTable(uda_path,working_dir)
  tau_yield = material_dict['PEAKI1']*material_dict['FSLOPE']
  #tau_yield = material_dict['PEAKI1']
  bulk_mod = material_dict['B0']
  shear_mod = material_dict['G0']
    
  analytical_times,analytical_sigmas,epsils=defTable_to_J2Solution(def_times,Fs,bulk_mod,shear_mod,tau_yield,num_substeps=1000)
   
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
    analytical_Sxx.append(sigma[0][0])
    analytical_Syy.append(sigma[1][1])
    analytical_Szz.append(sigma[2][2])

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$') 
  plt.figure(1)
  plt.clf()
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)    
  param_text = material_dict['material string']
  plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')    
  plt.title('TabularTest 11:\nUniaxial Strain J2 Plasticity (plot a)')
  plt.plot(np.array(analytical_e11),np.array(analytical_Sxx),'--g',linewidth=lineWidth+1,label=str_to_mathbf('Analytical'))    
  plt.plot(np.array(exx),np.array(Sxx),'-r',label=str_to_mathbf('Uintah'))
  plt.xlabel(str_to_mathbf('\epsilon_{A}'))
  plt.ylabel(str_to_mathbf('\sigma_{A} (Pa)'))     
  plt.legend()    
  savePNG(save_path+'/Test11_verificationPlot_a','1280x960')      
    
  plt.figure(2)
  plt.clf()
  ax2 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  ax2.xaxis.set_major_formatter(formatter)
  ax2.yaxis.set_major_formatter(formatter)    
  param_text = material_dict['material string']
  plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')     
  plt.title('TabularTest 11:\nUniaxial Strain J2 Plasticity (plot b)')    
  plt.plot(np.array(analytical_e11),np.array(analytical_Syy),'--g',linewidth=lineWidth+1,label=str_to_mathbf('Analytical'))
  plt.plot(np.array(exx),np.array(Syy),'-r',label=str_to_mathbf('Uintah'))
  plt.xlabel(str_to_mathbf('\epsilon_{A}'))
  plt.ylabel(str_to_mathbf('\sigma_{L} (Pa)')) 
  plt.legend()
  savePNG(save_path+'/Test11_verificationPlot_b','1280x960')      

  plt.figure(3)
  plt.clf()
  ax3 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  ax3.xaxis.set_major_formatter(formatter)
  ax3.yaxis.set_major_formatter(formatter)    
  param_text = material_dict['material string']
  plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')   
  plt.title('TabularTest 11:\nUniaxial Strain J2 Plasticity (plot c)')     
  plt.plot(analytical_times,np.array(analytical_e11),'-g',linewidth=lineWidth+1,label=str_to_mathbf('Analytical \epsilon_{xx}'))
  plt.plot(analytical_times,np.array(analytical_e22),'-r',linewidth=lineWidth+1,label=str_to_mathbf('Analytical \epsilon_{yy}'))
  plt.plot(analytical_times,np.array(analytical_e33),'-b',linewidth=lineWidth+1,label=str_to_mathbf('Analytical \epsilon_{zz}'))
  plt.legend()
  plt.xlabel(str_to_mathbf('Time (s)'))
  plt.ylabel(str_to_mathbf('\epsilon'))
  savePNG(save_path+'/Test11_verificationPlot_c','1280x960') 
    
  plt.figure(4)
  plt.clf()
  ax4 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  ax4.xaxis.set_major_formatter(formatter)
  ax4.yaxis.set_major_formatter(formatter)    
  param_text = material_dict['material string']
  plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')   
  plt.title('TabularTest 11:\nUniaxial Strain J2 Plasticity (plot d)')     
  plt.plot(analytical_times,np.array(analytical_Sxx),'-g',linewidth=lineWidth+1,label=str_to_mathbf('Analytical \sigma_{xx}'))
  plt.plot(analytical_times,np.array(analytical_Syy),'-r',linewidth=lineWidth+1,label=str_to_mathbf('Analytical \sigma_{yy}'))
  plt.plot(analytical_times,np.array(analytical_Szz),'-b',linewidth=lineWidth+1,label=str_to_mathbf('Analytical \sigma_{zz}'))
  plt.legend()
  plt.xlabel(str_to_mathbf('Time (s)'))
  plt.ylabel(str_to_mathbf('\sigma (Pa)'))    
  savePNG(save_path+'/Test11_verificationPlot_d','1280x960') 
    
  if SHOW_ON_MAKE:
    plt.show()
  

def test12_postProc(uda_path,save_path,**kwargs):
  #Extract stress history
  print("Post Processing Test: 12 - Nonlinear Elasticity")
  times,sigmas = get_pStress(uda_path)
  pressure = []
  for sigma in sigmas:
    pressure.append(-sigma_I1(sigma)/3.0)
  times,plasticStrainVol = get_pPlasticStrainVol(uda_path)
  times,elasticStrainVol = get_pElasticStrainVol(uda_path)
  totalStrainVol = -np.array(elasticStrainVol)-np.array(plasticStrainVol)

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
  plt.plot(totalStrainVol,pressure,'-b',label='Tabular')
  plt.title('TabularTest 12:\nNonlinear Elasticity')
  plt.ylabel(str_to_mathbf('p: pressure (Pa)'))
  plt.xlabel(str_to_mathbf('ev: compressive volumetric strain'))
  #ax1.set_xticks([0,0.005,0.010,0.015,0.020,0.025])
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)   
  plt.legend()
  savePNG(save_path+'/Test12_verificationPlot','1280x960')
  if SHOW_ON_MAKE:
    plt.show()

def test13_postProc(uda_path,save_path,**kwargs):
  COLORS = ['Black','Blue','Magenta','Red','Green']
  if 'WORKING_PATH' in kwargs:
    working_dir = kwargs['WORKING_PATH']  
  else:
    print('\nERROR: need working directory to post process this problem')
    return

  #Plot Constants
  Xlims = (-450,50)
  Ylims = (-100,100)  
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')  
  plt.figure(1)
  plt.hold(True)
  plt.clf()

  material_dict = get_yield_surface_data(uda_path)    
  PEAKI1 = material_dict['PEAKI1']
  FSLOPE = material_dict['FSLOPE']
  STREN = material_dict['STREN']
  T1 = material_dict['T1']
  T2 = material_dict['T2']    

  def_times,Fs = get_defTable(uda_path,working_dir)
  A = Fs[1][0][0]
  #As = Fs[10][0][0]
  K = material_dict['B0']
  G = material_dict['G0']
  C = K+(4.0/3.0)*G
  Y = STREN*1.732
  YS = STREN
    
  #uniaxial strain (unscaled)
  analytical_exx = [0.0, (Y/(2.0*G)), np.log(A), ]  
  analytical_Sxx=[0.0,
    (C*Y)/(2.0*G),
    (((C-K)*Y)/(2*G)+K*np.log(A)),
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
  times,sigmas = get_pStress(uda_path)
  times,epsils = get_epsilons(uda_path)
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
    Sxx.append(sigma[0][0])
    Syy.append(sigma[1][1])
    Szz.append(sigma[2][2])
    Sxy.append(sigma[0][1])
      
  scaled_exx = ((2.0*G)/Y)*np.array(exx)
  scaled_Sxx = ((2.0*G)/(C*Y))*np.array(Sxx)*1.0e6
  scaled_Syy = ((2.0*G)/(C*Y))*np.array(Syy)*1.0e6
  #S = np.array(Sxx) - np.array(Syy)
  S = np.array(Sxx)
  #E = np.array(exy)

  ###PLOTTING
  ax1 = plt.subplot(111)
  plt.subplots_adjust(right=0.75)
  #param_text = material_dict['material string']
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='x-small')   
  eqShear_vs_meanStress(exx,S,LINE_LABEL = 'T1='+format(T1,'1.3e')+' T2='+format(T2,'1.3e'))
  #eqShear_vs_meanStress(E,S,LINE_LABEL = 'T1='+format(T1,'1.3e')+' T2='+format(T2,'1.3e'),COLOR=COLORS[idx])
  plt.plot(analytical_exx,analytical_Sxx,'--',color='Red',label='Analytical solution for rate independent case.')
  plt.title('TabularTest 13:')
  plt.ylabel(str_to_mathbf('\sigma_{xx}'))
  plt.xlabel(str_to_mathbf('\epsilon_{xx}'))  
  #plt.ylabel(str_to_mathbf('\sigma_{xy}'))
  #plt.xlabel(str_to_mathbf('\epsilon_{xy}'))
  ax1.xaxis.set_major_formatter(formatter)
  ax1.yaxis.set_major_formatter(formatter)   
  plt.legend()
  savePNG(save_path+'/Test13_verificationPlot','1280x960')
  if SHOW_ON_MAKE:
    plt.show()   
  
#-----------------------------------------------------------------------------------
# Read the simulation stress data and compute p,q (converts to Pa)
#-----------------------------------------------------------------------------------
def readSimStressData(uda_path):

  NAN_FAIL = False

  #Extract stress history
  print("Extracting stress history...")
  args = [partextract_exe, "-partvar","p.stress",uda_path]
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

    sigma_a = -S11
    sigma_r = -S22
    sigma_ar = -S12
    
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
  for sigma in sigma_sim:
    I1_sim.append(sigma_I1(sigma))
    J2_sim.append(sigma_J2(sigma))

  # Compute p, q
  pp_sim = -np.array(I1_sim)/3.0
  qq_sim = np.sqrt(3.0*np.array(J2_sim))

  return time_sim, sigma_sim, sigma_a_sim, sigma_r_sim, sigma_ar_sim, pp_sim, qq_sim

#---------------------------------------------------------------------------------
# Interpolate the data at a set of given time points
#---------------------------------------------------------------------------------
def getDataTimeSnapshots(time_snapshots, time, data):

  # Create a clean list containing the data as functions of time
  time_list = []
  val_list = []
  for ii in range(1, len(time)):
    tt_0 = time[ii-1]
    tt_1 = time[ii]
    for jj, ta in enumerate(time_snapshots):
      ss = (ta - tt_0)/(tt_1 - tt_0)
      if (ss >= 0.0 and ss < 1.0):
        #print("ii = " , ii, " tt_0 = " , tt_0, " tt_1 = ", tt_1, " jj = " , jj, " ta = " , ta )
        time_list.append(ta)
        val_list.append((1-ss)*data[ii-1]+ss*data[ii]) 
   
  return time_list, val_list

#-----------------------------------------------------------------------------------
# Plot the sim data as a function of time (assume Pa)
#-----------------------------------------------------------------------------------
def plotSimDataSigmaTime(fig, time_snapshots, time_sim, sigma_a_sim, sigma_r_sim, sigma_ar_sim):

  # Get snapshots from data
  time_snap, sigma_a_snap = getDataTimeSnapshots(time_snapshots, time_sim, sigma_a_sim)
  time_snap, sigma_r_snap = getDataTimeSnapshots(time_snapshots, time_sim, sigma_r_sim)
  time_snap, sigma_ar_snap = getDataTimeSnapshots(time_snapshots, time_sim, sigma_ar_sim)

  # Activate the figure
  plt.figure(fig.number)

  # Plot sigma_a vs. time
  plt.plot(time_sim, sigma_a_sim, '--r', label='$\sigma_a$ (sim)')
  plt.plot(time_sim, sigma_r_sim, '--b', label='$\sigma_r$ (sim)')
  plt.plot(time_sim, sigma_ar_sim, '--g', label='$\sigma_{ar}$ (sim)')

  # Plot filled circles at time snapshots
  for ii in range(0, len(time_snap)):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii)/len(time_snap))
    plt.plot(time_snap[ii], sigma_a_snap[ii], 'o', color=plt_color) 
    plt.plot(time_snap[ii], sigma_r_snap[ii], 'o', color=plt_color) 
    #plt.plot(time_snap[ii], sigma_ar_snap[ii], 'o', color=plt_color) 

  return time_snap, sigma_a_snap, sigma_r_snap, sigma_ar_snap

#-----------------------------------------------------------------------------------
# Plot the sim data (pq) as a function of time (aasume Pa)
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
    plt_color = cm.Paired(float(ii)/len(time_snap))
    plt.plot(time_snap[ii], p_snap[ii], 'o', color=plt_color) 
    plt.plot(time_snap[ii], q_snap[ii], 'o', color=plt_color) 

  return time_snap, p_snap, q_snap

