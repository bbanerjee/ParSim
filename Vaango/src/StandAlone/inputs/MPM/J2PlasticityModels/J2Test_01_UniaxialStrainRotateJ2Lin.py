from J2TestSuite_PostProcUtils import *
from J2YieldSurfaceUtils import *


def uniaxialStrainRotateJ2Lin(uda_path, save_path, **kwargs):
  print("Post Processing Test: 01 - Uniaxial Compression With Rotation")

  # Read the stress simulation data
  times, sigmas, sigma_a_sim, sigma_r_sim, sigma_ar_sim, pp_sim, qq_sim = readSimStressData(
      uda_path)

  # Set up rotation angles and matrices
  rotation_angles = list(
      map(lambda t, tmax: t * 360 / tmax, times, [max(times)] * len(times)))
  cosines = list(map(lambda t: np.cos(t * math.pi / 180), rotation_angles))
  sines = list(map(lambda t: np.sin(t * math.pi / 180), rotation_angles))
  rot_mats = list(
      map(lambda c, s: np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]]), cosines,
          sines))
  #print(rot_mats)

  # Compute the stress components in the rotated coordinate system
  sigmas_rot = list(
      map(lambda Q, sigma: np.dot(np.dot(Q, sigma), Q.transpose()), rot_mats,
          sigmas))
  #print(sigmas_rot)

  # Set up time points
  analytical_times = np.linspace(0.0, times[-1], 15)

  # Read the simulation state data
  ep_sim, epdot_sim, backstress_list, phi_list, D_list, T_sim, vol_sim, time_list = getInternalVariables(
      uda_path, analytical_times)

  # Read the deformation data
  times, defGrad_sim = get_pTensor(uda_path, "p.deformationGradinet", 0)
  J_sim = [tensor_det(F) for F in defGrad_sim]

  # Get the model parameters
  material_dict = get_yield_surface_data(uda_path)
  param_text = material_dict['material string']

  # Get snapshots of pq data
  t_sim_snap, p_sim_snap = getDataTimeSnapshots(analytical_times, times, pp_sim)
  t_sim_snap, q_sim_snap = getDataTimeSnapshots(analytical_times, times, qq_sim)

  # Get snapshots of state data
  t_sim_snap, ep_sim_snap = getDataTimeSnapshots(analytical_times, times,
                                                 ep_sim)
  t_sim_snap, epdot_sim_snap = getDataTimeSnapshots(analytical_times, times,
                                                    epdot_sim)
  t_sim_snap, T_sim_snap = getDataTimeSnapshots(analytical_times, times, T_sim)
  t_sim_snap, J_sim_snap = getDataTimeSnapshots(analytical_times, times,
                                                  J_sim)
  rho_0 = material_dict['rho0']
  rho_sim_snap = [rho_0 * J for J in J_sim_snap]

  # Find the plot limits
  Sxx_rot = [sigma[0][0] for sigma in sigmas_rot]
  Syy_rot = [sigma[1][1] for sigma in sigmas_rot]
  Sxy_rot = [sigma[0][1] for sigma in sigmas_rot]

  ###PLOTTING
  formatter = ticker.FormatStrFormatter('$\mathbf{%g}$')
  param_text = material_dict['material string']

  #---------------------------------------------------------------------------------
  # Plot rotated stress as a function of time
  fig1 = plt.figure(1)
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='xx-small')
  plotSimDataSigmaTime(fig1, analytical_times, times, Sxx_rot, Syy_rot, Sxy_rot,
                       '$\sigma_{a}$ (sim)', '$\sigma_{r}$ (sim)',
                       '$\sigma_{ar}$ (sim)')

  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  plt.xlabel(str_to_mathbf('Time (sec)'))
  plt.ylabel(str_to_mathbf('Rotated Stress (Pa)'))
  plt.grid(True)
  plt.legend(loc='best', prop={'size': 10})
  savePNG(save_path + '/UniaxialStrainRotateJ2Lin_sigma_rot_time', '1280x960')

  #----------------------------------------------------------------
  # Plot the yield surface for test1
  #----------------------------------------------------------------
  # Set up figure
  fig2 = plt.figure(2)
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='xx-small')

  # Plot p vs. q simulation results
  eqShear_vs_meanStress(pp_sim, qq_sim)

  # Plot filled circles at time snapshots
  for ii in range(0, len(t_sim_snap) - 1):

    # Choose the Paired colormap
    plt_color = cm.Paired(float(ii) / len(t_sim_snap))
    plt.plot(p_sim_snap[ii], q_sim_snap[ii], 'o', color=plt_color)

    # Plot yield surfaces
    plotPQYieldSurfaceSim(plt, material_dict, pp_sim, epdot_sim_snap[ii],
                          ep_sim_snap[ii], T_sim_snap[ii], rho_sim_snap[ii])

  savePNG(save_path + '/UniaxialStrainRotateJ2Lin_yield_surface', '1280x960')
  #plt.show()

  #---------------------------------------------------------------------------------
  # Plot experimental and simulation data as a function of time
  fig3 = plt.figure(3)
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='xx-small')
  plotSimDataSigmaTime(fig3, analytical_times, times, sigma_a_sim, sigma_r_sim,
                       sigma_ar_sim, '$\sigma_{xx}$ (sim)',
                       '$\sigma_{yy}$ (sim)', '$\sigma_{xy}$ (sim)')

  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  plt.xlabel(str_to_mathbf('Time (sec)'))
  plt.ylabel(str_to_mathbf('Stress (Pa)'))
  plt.grid(True)
  plt.legend(loc='best', prop={'size': 10})
  savePNG(save_path + '/UniaxialStrainRotateJ2Lin_sigma_time', '1280x960')

  fig4 = plt.figure(4)
  plt.clf()
  #plt.subplots_adjust(right=0.75)
  #plt.figtext(0.77,0.70,param_text,ha='left',va='top',size='xx-small')
  plotSimDataPQTime(fig4, analytical_times, times, pp_sim, qq_sim)
  axes = plt.gca()
  axes.xaxis.set_major_formatter(formatter)
  axes.yaxis.set_major_formatter(formatter)
  plt.xlabel(str_to_mathbf('Time (sec)'))
  plt.ylabel(str_to_mathbf('Stress (Pa)'))
  plt.grid(True)
  plt.legend(loc='best', prop={'size': 8})
  savePNG(save_path + '/UniaxialStrainRotateJ2Lin_pq_time', '1280x960')

  plt.show()
