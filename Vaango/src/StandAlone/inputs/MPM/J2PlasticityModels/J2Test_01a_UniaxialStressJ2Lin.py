from TabularTestSuite_PostProcUtils import *
from TabularYieldSurfaceUtils import *

def uniaxialStressJ2Lin(uda_path, save_path,**kwargs):
  print("Post Processing Test: 01a - Uniaxial Stress Compression")

  # Read the stress simulation data
  times, sigmas, sigma_a_sim, sigma_r_sim, sigma_ar_sim, pp_sim, qq_sim = readSimStressData(uda_path, matID = 1)

  # Set up time points
  analytical_times = np.linspace(0.0, times[-1], 15)

  # Read the interval variable simulation data
  ev_e_list, ev_p_list, times_list, ev_e, ev_p = getInternalVariables(uda_path, analytical_times, matID = 1)

  # Get the model parameters
  material_dict = get_yield_surface_data(uda_path)
  param_text = material_dict['material string']
  elastic_table = getJSONTable(material_dict['elastic_filename'])
  yield_table = getJSONTable(material_dict['yield_filename'])

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
  print("Sxx_min = ", Sxx_min)
  print("Sxx_max = ", Sxx_max)
  print("Syy_min = ", Syy_min)
  print("Syy_max = ", Syy_max)

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$') 
  param_text = material_dict['material string']

  #----------------------------------------------------------------
  # Plot the yield surface for test1
  #----------------------------------------------------------------
  # Set up figure
  fig1 = plt.figure(1)
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='xx-small')  

  # Plot p vs. q simulation results
  eqShear_vs_meanStress(pp_sim, qq_sim)  

  # Plot filled circles at time snapshots
  for ii in range(0, len(t_sim_snap)-1):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii)/len(t_sim_snap))
    plt.plot(p_sim_snap[ii], q_sim_snap[ii], 'o', color=plt_color) 

  # Plot yield surfaces
  pmin = min(pp_sim)
  pmax = max(pp_sim)
  qmax = max(map(lambda q : abs(q), qq_sim))
  plotPQYieldSurfaceSim(plt, material_dict, yield_table,
                        ev_e_list, ev_p_list, times_list,
                        pmin, pmax, qmax) 

  savePNG(save_path+'/UniaxialStressJ2Lin_yield_surface','1280x960')
  #plt.show()

  #---------------------------------------------------------------------------------
  # Plot experimental and simulation data as a function of time
  fig2 = plt.figure(2)
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='xx-small')  
  plotSimDataSigmaTime(fig2, analytical_times, times, sigma_a_sim, sigma_r_sim, sigma_ar_sim,
                       '$\sigma_{xx}$ (sim)', '$\sigma_{yy}$ (sim)',
                       '$\sigma_{xy}$ (sim)')

  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  plt.xlabel(str_to_mathbf('Time (sec)')) 
  plt.ylabel(str_to_mathbf('Stress (Pa)')) 
  plt.grid(True)
  plt.legend(loc='best', prop={'size':10}) 
  savePNG(save_path+'/UniaxialStressJ2Lin_sigma_time','1280x960')

  fig3 = plt.figure(3)
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='xx-small')  
  plotSimDataPQTime(fig3, analytical_times, times, pp_sim, qq_sim)
  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  plt.xlabel(str_to_mathbf('Time (sec)')) 
  plt.ylabel(str_to_mathbf('Stress (Pa)')) 
  plt.grid(True)
  plt.legend(loc='best', prop={'size':8}) 
  savePNG(save_path+'/UniaxialStressJ2Lin_pq_time','1280x960')

  plt.show()



