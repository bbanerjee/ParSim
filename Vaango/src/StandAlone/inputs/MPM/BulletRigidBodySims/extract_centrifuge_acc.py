#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import multiprocessing
import subprocess as sub_proc  
import tempfile
import argparse
import numpy as np

#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib import ticker

#--------------------------------------------------------------------------
# Main function (parse arguments)
#--------------------------------------------------------------------------
def main(argv):

  uda_name = ''
  material_id = 0 

  # Parse comman line arguments
  parser = argparse.ArgumentParser(
    usage="extract_centrifuge_acc [--x|--y] --point <x> <y> <z> --mat <mat id> --uda <uda file_name_prefix>")
  parser.add_argument('--x', action='store_const', dest='direction', const='x', required=False)
  parser.add_argument('--y', action='store_const', dest='direction', const='y', required=False)
  parser.add_argument('--point', nargs=3, help = 'x y z', type=float, required=True)  # For z-value
  parser.add_argument('--mat', default=0, type=int, required=True)
  parser.add_argument('--uda', required=True)
  args = parser.parse_args()

  vals = args.point
  if (vals):
    x = float(vals[0])
    y = float(vals[1])
    z = float(vals[2])

  # Parse box
  if (args.direction == 'x'):
    x = -99999
    #y = 0

  if (args.direction == 'y'):
    #x = 0
    y = -99999


  # Parse mat
  material_id = int(args.mat)

  # Parse uda
  uda_name = args.uda
  
  return uda_name, x, y, z, material_id

#--------------------------------------------------------------------------
# Define the executables to be used and the working directories 
# Two executables are needed - selectpart and partextract
# selectpart:  selects a list of particleIDs inside a selection box
# partextract: extracts particle variables corresponding to each particleID 
#              in the selection box
# In bash:
#   export SELECTPART_EXE=/home/..../selectpart
#   export PARTEXTRACT_EXE=/home/..../partextract
#--------------------------------------------------------------------------
def set_execs_and_dirs():

  # Define the executables
  if (os.getenv('TIMEXTRACT_EXE') == None):
    print "**ERROR** Missing environment var. In bash: export TIMEXTRACT_EXE=/path/to/timeextract"
    sys.exit()
  extract_exe = os.path.abspath(os.environ['TIMEXTRACT_EXE'])

  # Define the working and output directories
  # The uda files will be read from the working directory
  if (os.getenv('WORKING_DIR') == None):
    print "**ERROR** Missing environment var. In bash: export WORKING_DIR=/path/to/uda/directory"
    sys.exit()
  work_dir = os.path.abspath(os.environ['WORKING_DIR'])

  # Make a directory for saving the plots and extracted data
  plot_dir = os.path.join(work_dir, "DataAndPlots")
  if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

  return extract_exe, work_dir, plot_dir


#--------------------------------------------------------------------------
# Run timeextract for each particle id in the selected list
#--------------------------------------------------------------------------
def timeextract(extract_exe, work_dir, plot_dir, uda_name, x, y, z, material_id):

  # Create the full uda path
  uda_dir = os.path.join(work_dir, uda_name)
  if not os.path.exists(uda_dir):
    print "**ERROR** Uda file does not exist", uda_dir
    sys.exit()

  # Create the command + arguments (for the position)
  args_pos = [extract_exe, "-p", str(x), str(y), str(z), 
    "-v", "g.acceleration", "-uda", uda_dir]
  print args_pos

  # Run the command and save output in part_data_files
  pos_data = tempfile.TemporaryFile()
  error_file = os.path.join(plot_dir, uda_name + ".error")
  err_file = open(error_file+".time_extract", "w+")
  #tmp = sub_proc.Popen(args_pos, stdout=pos_data, stderr=sub_proc.PIPE)
  tmp = sub_proc.Popen(args_pos, stdout=pos_data, stderr=err_file)
  dummy = tmp.wait()

  # Read the part data file to get a list of positions as a function of time
  pos_data.seek(0)
  times = []
  ax = []
  ay = []
  az = []
  for line in pos_data:
    line = line.replace("[", "").replace("]","") 
    line = line.strip().split()
    print line
    try:
      float(line[0])
    except ValueError:
      continue
    times.append(float(line[0]))
    ax.append(np.float64(line[1]))
    ay.append(np.float64(line[2]))
    az.append(np.float64(line[3]))

  # close and return
  pos_data.close()
  err_file.close()

  # Join all the data into one data frame with columns
  # time, ax, ay, az
  acc_data = zip(times, ax, ay, az)

  return acc_data

#--------------------------------------------------------------------------
# Set up plot parameters
#--------------------------------------------------------------------------
def set_plot_parameters():

  #Set matplotlib defaults to desired values
  #S#Set the legend to best fit
  fontSize = 16
  markers = None
  plt.rcParams['legend.loc']='right'
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

#--------------------------------------------------------------------------
# Save as pdf
#--------------------------------------------------------------------------
def savePDF(fig, name, size='1920x1080'):
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
  fig.set_size_inches(size[0],size[1])
  #save at speciified resolution
  fig.savefig(name+'.pdf', bbox_inches=0, dpi=plt.rcParams['figure.dpi'])

def savePNG(fig, name, size='1920x1080'):
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
  fig.set_size_inches(size[0],size[1])
  #save at speciified resolution
  fig.savefig(name+'.png', bbox_inches=0, dpi=plt.rcParams['figure.dpi'])

#--------------------------------------------------------------------------
# For bold math fonts
#--------------------------------------------------------------------------
def str_to_mathbf(string):
  #Only works with single spaces no leading space
  string = string.split()
  return_string = ''
  for elem in string:
    elem = r'$\mathbf{'+elem+'}$'
    return_string+=elem+'  '
  return return_string[0:-1]

#--------------------------------------------------------------------------
# For exponential format
#--------------------------------------------------------------------------
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

#--------------------------------------------------------------------------
# Plot the acc as a function of time
#--------------------------------------------------------------------------
def plot_acc(fig, acc_data, direction, location, plt_color, line_style):

  # Function that returns t
  # Get the data into plottable form
  times = []
  acc = []
  for data in acc_data:
    times.append(data[0]+0.006)
    #acc.append(np.sqrt(data[1]*data[1] + data[2]*data[2] + data[3]*data[3])/9.81)
    if (direction == 'x'):
      acc.append(data[1]/9.81)
    elif (direction == 'y'):
      acc.append(data[2]/9.81)
    else:
      acc.append(data[3]/9.81)

  # Plot x-t
  plt.figure(fig)
  ax = plt.subplot(111)
  plt.plot(np.array(times), np.array(acc), linestyle=line_style,
    linewidth=2, color=plt_color, label='r = '+str(location*100)+' cm')

  plt.xlabel(str_to_mathbf('Time (s)'))
  plt.ylabel(str_to_mathbf(direction+'-Acceleration (g)'))

  formatter_int = ticker.FormatStrFormatter('$\mathbf{%g}$')
  formatter_exp = ticker.FuncFormatter(exp_fmt)

  #ax.xaxis.set_major_formatter(formatter_exp)
  #ax.yaxis.set_major_formatter(formatter_exp)    
  ax.xaxis.set_major_formatter(formatter_int)
  ax.yaxis.set_major_formatter(formatter_int)

  ax.set_xlim(0.005, 0.01)
  #ax.set_ylim(Ylims[0],Ylims[1])
  ax.grid(True)
  
  return ax

#-----------------------------------------------------------------------
# Acceleraometer positions
#-----------------------------------------------------------------------
def getPositions(direction, zpos):

  rad = [0.152, 0.203, 0.254, 0.305, 0.356]
  x = []
  y = []
  z = []
  for pos in rad:
    if (direction == 'x'):
      x.append(np.float64(pos))
      y.append(np.float64(0.0))
      z.append(np.float64(zpos))
    else:
      x.append(np.float64(0.0))
      y.append(np.float64(pos))
      z.append(np.float64(zpos))
     
  return zip(x, y, z)

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Main entry point
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
if __name__ == "__main__":

  uda_name, x, y, z, material_id = main(sys.argv[1:])
  print 'Input uda file is : ', uda_name
  print '  Point : ', x, y, z
  print '  Material id : ', material_id

  extract_exe, work_dir, plot_dir = set_execs_and_dirs()
  print extract_exe, work_dir, plot_dir

  # Set up figures
  set_plot_parameters()
  fig1 = plt.figure("x-acc")
  plt.clf()
  plt.subplots_adjust(right=0.75)
  fig2 = plt.figure("y-acc")
  plt.clf()
  plt.subplots_adjust(right=0.75)
  fig3 = plt.figure("z-acc")
  plt.clf()
  plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70, param_text,ha='left',va='top',size='x-small')

  if (x == -99999):
    positions = getPositions('x', z)
  elif (y == -99999):
    positions = getPositions('y', z)
  else:
    positions = [[x, y, z]]
  print positions

  colors = ['Red', 'Green', 'Blue', 'Magenta', 'Violet']
  for idx, pos in enumerate(positions):
    x = pos[0]
    y = pos[1]
    z = pos[2]
    rad = np.sqrt(x*x + y*y)
    data = timeextract(extract_exe, work_dir, plot_dir, uda_name, x, y, z, material_id)
    plt_color = colors[idx]
    line_style = '-'
    plot_acc('x-acc', data, 'x', rad, plt_color, line_style)
    plot_acc('y-acc', data, 'y', rad, plt_color, line_style)
    plot_acc('z-acc', data, 'z', rad, plt_color, line_style)

  plt.figure('x-acc')
  plt.legend(bbox_to_anchor=(1.05,1), loc=2, prop={'size':10})
  plt.figure('y-acc')
  plt.legend(bbox_to_anchor=(1.05,1), loc=2, prop={'size':10})
  plt.figure('z-acc')
  plt.legend(bbox_to_anchor=(1.05,1), loc=2, prop={'size':10})
  plt_file_x = os.path.join(plot_dir, uda_name+".xacc")
  plt_file_y = os.path.join(plot_dir, uda_name+".yacc")
  plt_file_z = os.path.join(plot_dir, uda_name+".zacc")
  savePNG(fig1, plt_file_x, size='1280x960')
  savePNG(fig2, plt_file_y, size='1280x960')
  savePNG(fig3, plt_file_z, size='1280x960')
  plt.show()


