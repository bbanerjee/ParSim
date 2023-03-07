from J2TestSuite_PostProcUtils import *
from J2YieldSurfaceUtils import *

def triaxialStrainTension(uda_path, save_path,**kwargs):
  print("Post Processing Test: 11 - Triaxial strain tension")

  # Read the stress simulation data
  times, sigmas, sigma_a_sim, sigma_r_sim, sigma_ar_sim, pp_sim, qq_sim = readSimStressData(uda_path, matID = 0)

  # Set up time points
  analytical_times = np.linspace(0.0, times[-1], 15)

  # Read the interval variable simulation data
  ev_e_list, ev_p_list, times_list, ev_e_sim, ev_p_sim = \
    getInternalVariables(uda_path, analytical_times, matID = 0)

  # Get the model parameters
  material_dict = get_yield_surface_data(uda_path)
  param_text = material_dict['material string']
  elastic_table = getJSONTable(material_dict['elastic_filename'])
  yield_table = getJSONTable(material_dict['yield_filename'])
  hydrostat_table = getJSONTable('DrySand_HydrostatData.json')

  # Extract the data from the hydrostat table
  ev_hydrostat   = hydrostat_table['TotalStrainVol']
  pbar_hydrostat = hydrostat_table['Pressure']

  # Get snapshots of pq data 
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
  #print("Sxx_min = ", Sxx_min)
  #print("Sxx_max = ", Sxx_max)
  #print("Syy_min = ", Syy_min)
  #print("Syy_max = ", Syy_max)

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$') 
  param_text = material_dict['material string']
  compression = 'positive'

  #----------------------------------------------------------------
  # Plot the yield surface 
  #----------------------------------------------------------------
  # Set up figure
  fig2 = plt.figure(2)
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='xx-small')  

  # Plot p vs. q simulation results
  eqShear_vs_meanStress(pp_sim, qq_sim, compression)

  # Plot filled circles at time snapshots
  for ii in range(0, len(t_sim_snap)):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii)/len(t_sim_snap))
    if (compression == 'positive'):
      plt.plot(-p_sim_snap[ii], -q_sim_snap[ii], 'o', color=plt_color) 
    else:
      plt.plot(p_sim_snap[ii], q_sim_snap[ii], 'o', color=plt_color) 

  # Plot yield surface
  pbarmin = min(yield_table['Pressure'])
  pbarmax = max(map(lambda p: abs(p), pp_sim))
  qmax = max(map(lambda q: abs(q), qq_sim))
  plotPQYieldSurfaceSim(plt, material_dict, yield_table,
                          ev_e_list, ev_p_list, times_list,
                          pbarmin, pbarmax, qmax, compression)

  savePNG(save_path+'/TriaxialStrainTension_yield_surface','1280x960')
  #plt.show()

  #---------------------------------------------------------------------------------
  # Plot experimental and simulation data as a function of time
  fig3 = plt.figure(3)
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='xx-small')  
  plotSimDataSigmaTime(fig3, analytical_times, times, sigma_a_sim, sigma_r_sim, sigma_ar_sim,
                       '$\sigma_{xx}$ (sim)', '$\sigma_{yy}$ (sim)',
                       '$\sigma_{xy}$ (sim)', compression)

  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  plt.xlabel(str_to_mathbf('Time (sec)')) 
  plt.ylabel(str_to_mathbf('Stress (Pa)')) 
  plt.grid(True)
  plt.legend(loc='best', prop={'size':10}) 
  savePNG(save_path+'/TriaxialStrainTension_sigma_time','1280x960')

  fig3 = plt.figure(3)
  plt.clf()

  plotSimDataSigmaEps(fig3, analytical_times, times, pp_sim, ev_e_sim, ev_p_sim,
                      compression = 'positive')
  plt_color = cm.Paired(1)
  #plt.plot(ev_hydrostat, pbar_hydrostat, '-', color=plt_color, label='Experimental data') 

  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  #axes.set_xlim([0, 0.5])
  #axes.set_ylim([0, 1.2*max(pbar_hydrostat)])
  plt.xlabel(str_to_mathbf('Strain ')) 
  plt.ylabel(str_to_mathbf('Stress (Pa)')) 
  plt.grid(True)
  plt.legend(loc='best', prop={'size':10}) 
  savePNG(save_path+'/TriaxialStrainTension_pbar_evbar','1280x960')
  #plt.show()

  fig4 = plt.figure(4)
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='xx-small')  
  plotSimDataPQTime(fig4, analytical_times, times, pp_sim, qq_sim, compression)
  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  plt.xlabel(str_to_mathbf('Time (sec)')) 
  plt.ylabel(str_to_mathbf('Stress (Pa)')) 
  plt.grid(True)
  plt.legend(loc='best', prop={'size':8}) 
  savePNG(save_path+'/TriaxialStrainTension_pq_time','1280x960')

  plt.show()



