import numpy as np
import pyvista as pv

def principal_stresses(stress):

  I = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).reshape(3,3)
  I1 = stress.trace()
  s = stress - I * (I1 / 3.0)
  ss = np.matmul(s, s)
  J2 = ss.trace() / 2.0
  J3 = np.linalg.det(s)
  p = I1 / 3.0
  q = np.sqrt(3.0 * J2)
  r_cubed = 27.0/2.0 * J3
  cos3theta = r_cubed / (q * q * q)
  if (cos3theta < -1):
    cos3theta = -1
  if (cos3theta > 1):
    cos3theta = 1
  theta = np.arccos(cos3theta) / 3.0

  sigma_1 = p + 2.0 / 3.0 * q * np.cos(theta)
  sigma_2 = p + 2.0 / 3.0 * q * np.cos(theta - 2.0 * np.pi / 3.0)
  sigma_3 = p + 2.0 / 3.0 * q * np.cos(theta + 2.0 * np.pi / 3.0)

  sig_prin = np.array([sigma_1, sigma_2, sigma_3])

  return sig_prin

def lode_coordinates(stress):

  I = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).reshape(3,3)
  I1 = stress.trace()
  s = stress - I * (I1 / 3.0)
  ss = np.matmul(s, s)
  J2 = ss.trace() / 2.0
  J3 = np.linalg.det(s)
  p = I1 / 3.0
  q = np.sqrt(3.0 * J2)
  z = I1 / np.sqrt(3)
  r = np.sqrt(2.0 * J2)
  r_cubed = 27.0/2.0 * J3
  cos3theta = r_cubed / (q * q * q)
  if (cos3theta < -1):
    cos3theta = -1
  if (cos3theta > 1):
    cos3theta = 1
  theta = np.arccos(cos3theta) / 3.0

  #print("theta = ", theta * 180 / np.pi, "cos3theta = ", cos3theta, "cos(3theta) = ", np.cos(3.0*theta), "J3 = ", J3)
  #xx = z
  #yy = r * np.cos(theta) * np.sign(J3)
  #zz = r * np.sin(theta) * np.sign(J3)

  xx = z
  yy = r * np.cos(theta)
  zz = r * np.sin(theta)

  return np.array([xx, yy, zz])
  
def toMatrix(vec):
  return np.array([vec[0], vec[3], vec[4], vec[3], vec[1], vec[5], vec[4], vec[5], vec[2]]).reshape(3,3)

def rotateToInit(eigvec, stress):
  return np.matmul(np.matmul(eigvec.transpose(), stress), eigvec)

def I1J2(stress):
  I = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).reshape(3,3)
  I1 = stress.trace()
  s = stress - I * (I1 / 3.0)
  ss = np.matmul(s, s)
  J2 = ss.trace() / 2.0
  return I1, J2
  
def computePn(df_dsigma):

  k = 9.0e7
  mu = 6.0e7
  llambda = k - 2.0 / 3.0 * mu
  I = np.array([1.0, 1.0, 1.0, 0, 0, 0])
  tr_df_dsigma = df_dsigma[0] + df_dsigma[1] + df_dsigma[2]
  Pn = llambda * tr_df_dsigma * I + 2 * mu * df_dsigma
  Pn = Pn / np.linalg.norm(Pn)
  return Pn
  

def mohr_coulomb_lode(c, phi, low, high, numcells):

    rad = np.linspace(0.01, 1.73*high, numcells)
    grid = pv.CylinderStructured(radius = rad, height = np.sqrt(3)*(high - low),
                                 center = [0.5*(high+low), 0, 0],
                                 direction = [1, 0, 0],
                                 theta_resolution = 180,
                                 z_resolution = np.int(1.73*numcells))

    vtkGrid = pv.UnstructuredGrid(grid)
    stresses = vtkGrid.points

    sin_phi = np.sin(phi * np.pi / 180.0)
    cos_phi = np.cos(phi * np.pi / 180.0)
    sqrt3 = np.sqrt(3)

    yield_fn = []
    for stress in stresses:

      z = stress[0]
      r = np.sqrt(stress[1] * stress[1] + stress[2] * stress[2])
      Theta_bar = np.arctan2(stress[2], stress[1])
      theta = 1.0 / 3.0 * np.arccos(np.cos(3.0 * Theta_bar))
      #print("theta = ", theta * 180 / np.pi, "theta_bar = ", Theta_bar * 180 / np.pi)

      p = 1.0 / sqrt3 * z
      q = np.sqrt(3.0 / 2.0) * r

      Theta = theta + np.pi / 3.0
      Rmc = 1 / sqrt3 * np.sin(Theta) - 1 / 3 * np.cos(Theta) * sin_phi
      f = Rmc * q - p * sin_phi - c * cos_phi

      sigma_1 = p + 2.0 / 3.0 * q * np.cos(theta)
      sigma_2 = p + 2.0 / 3.0 * q * np.cos(theta - 2.0 * np.pi / 3.0)
      sigma_3 = p + 2.0 / 3.0 * q * np.cos(theta + 2.0 * np.pi / 3.0)

      ind = np.argsort(np.array([sigma_1, sigma_2, sigma_3]))
      i1 = ind[2]
      i3 = ind[0]
      f = f / (abs(stress[i1]) + abs(stress[i3]) + 2.0 * c)

      #print("z = ", z, "r = ", r, "theta = ", theta * 180 / np.pi, "f = ", f)
      yield_fn.append(f)

    vtkGrid["data"] = np.array(yield_fn)
    return vtkGrid

def drawArrow(start, end, scale):

  length = np.linalg.norm(end - start)
  direction = (end - start) / length
  cyl_bot = start
  cyl_top = end - direction * length * 0.3
  cyl_cen = 0.5 * (cyl_bot + cyl_top)
  arrow_cyl = pv.Cylinder(cyl_cen, direction = direction, radius = scale, 
                          height = 0.7 * length)
  cone_bot = cyl_top
  cone_top = cone_bot + direction * length * 0.3
  cone_cen = 0.5 * (cone_bot + cone_top)
  arrow_cone = pv.Cone(cone_cen, direction, height = np.linalg.norm(cone_top - cone_bot), 
                       radius = 2.5*scale)

  return arrow_cyl, arrow_cone
  
def calcData(p, df_dsigma, P_n, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad):

  df_dsigma = toMatrix(df_dsigma)
  P_n = toMatrix(P_n)

  alpha = np.linalg.norm(stress_trial - stress_alpha)

  stress_old = stress_old.reshape(3,3)
  stress_trial = stress_trial.reshape(3,3)
  stress_alpha = stress_alpha.reshape(3,3)

  stress_old_lode = lode_coordinates(stress_old)
  stress_trial_lode = lode_coordinates(stress_trial)
  stress_alpha_lode = lode_coordinates(stress_alpha)

  stress_old_ball = pv.Sphere(radius=ball_rad, center = stress_old_lode)
  stress_trial_ball = pv.Sphere(radius=ball_rad, center = stress_trial_lode)
  stress_alpha_ball = pv.Sphere(radius=ball_rad, center = stress_alpha_lode)

  df_dsigma_lode = lode_coordinates(df_dsigma)
  P_n_lode = lode_coordinates(P_n)
  stress_Pn_lode = stress_trial_lode - P_n_lode * arrow_len
  stress_df_lode = stress_old_lode + df_dsigma_lode * arrow_len
  #print("P_n_lode = ", P_n_lode)

  #stress_Pn = stress_trial - P_n * arrow_len 
  #stress_df = stress_old + df_dsigma * arrow_len
  #stress_Pn_lode = lode_coordinates(stress_Pn)
  #stress_df_lode = lode_coordinates(stress_df)

  stress_Pn_ball = pv.Sphere(radius=ball_rad, center = stress_Pn_lode)
  stress_df_ball = pv.Sphere(radius=ball_rad, center = stress_df_lode)

  df_dsigma_cyl, df_dsigma_cone = drawArrow(stress_old_lode, stress_df_lode, cyl_rad)
  P_n_cyl, P_n_cone = drawArrow(stress_trial_lode, stress_Pn_lode, cyl_rad)
  sigma_cyl, sigma_cone = drawArrow(stress_old_lode, stress_trial_lode, cyl_rad)

  blue = np.array([12/256, 238/256, 246/256])
  grey = np.array([189/256, 189/256, 189/256])
  yellow = np.array([255/256, 247/256, 0/256])

  p.add_mesh(stress_old_ball, color="red", show_edges=False)
  p.add_mesh(stress_trial_ball, color=blue, show_edges=False)
  p.add_mesh(stress_alpha_ball, color="#A4C639", show_edges=False)
  p.add_mesh(df_dsigma_cyl, color="green", show_edges=False)
  p.add_mesh(df_dsigma_cone, color="green", show_edges=False)
  p.add_mesh(P_n_cyl, color=yellow, show_edges=False)
  p.add_mesh(P_n_cone, color=yellow, show_edges=False)
  p.add_mesh(sigma_cyl, color="blue", show_edges=False)
  p.add_mesh(sigma_cone, color="blue", show_edges=False)

  return p

def plotData():

  fac = 0.1
  upper_lim = 1.0e5 * fac
  ball_rad = 2000 * fac
  cyl_rad = 1000 * fac
  arrow_len = 1.0e5 * fac 

  c = 1.0e4
  phi = 30
  grid = mohr_coulomb_lode(c, phi, -1.0e5, upper_lim, 30)
  contours = grid.contour([0.0])

  pv.set_plot_theme('document')
  p = pv.Plotter(notebook=False)
  #p.add_mesh(grid, opacity=0.3)
  p.add_mesh(contours, color="#FAEBD7", opacity=1.0, specular=1.0, show_scalar_bar=False)
  p.add_axes()
  #p.add_bounding_box()

  #df_dsigma = np.array([-0.0424937, -0.914076, 0.301887, 0.26744, 1.32117e-07, -9.68111e-08])
  #P_n = np.array([-0.250789, -0.944093, 0.0231495, 0.212736, 1.05093e-07, -7.70088e-08])
  #stress_old   = np.array([-9296.9, 1110.57, 0.00164243, 1110.57, -12916.2, -0.000385837, 0.00164243, -0.000385837, -5048.44])
  #stress_trial = np.array([-24293.1, 2901.95, 0.00429172, 2901.95, -33750.5, -0.00100821, 0.00429172, -0.00100821, -13191.7])
  #stress_alpha = np.array([-17524.2, -2839.85, 0.00145523, -2839.85, -8269.17, 0.00107028, 0.00145523, 0.00107028, -13816.5])
  #p = calcData(p, df_dsigma, P_n, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)
  
  #df_dsigma = np.array([-0.0424937, -0.914076, 0.301887, 0.26744, 1.32117e-07, -9.68111e-08])
  #P_n = np.array([-0.250789, -0.944093, 0.0231495, 0.212736, 1.05093e-07, -7.70088e-08])
  #stress_old   = np.array([-9296.9, 1110.57, 0.00164243, 1110.57, -12916.2, -0.000385837, 0.00164243, -0.000385837, -5048.44])
  #stress_trial = np.array([-16795, 2006.26, 0.00296707, 2006.26, -23333.4, -0.000697021, 0.00296707, -0.000697021, -9120.08])
  #stress_alpha = np.array([-13410.6, -864.644, 0.00154883, -864.644, -10592.7, 0.000342224, 0.00154883, 0.000342224, -9432.49])
  #p = calcData(p, df_dsigma, P_n, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.0424937, -0.914076, 0.301887, 0.26744, 1.32117e-07, -9.68111e-08])
  #P_n = np.array([-0.250789, -0.944093, 0.0231495, 0.212736, 1.05093e-07, -7.70088e-08])
  #stress_old   = np.array([-14486.9, 48.3381, 0.00199985, 48.3381, -14644.4, 1.1731e-05, 0.00199985, 1.1731e-05, -9333.14])
  #stress_trial = np.array([-21984.9, 944.028, 0.00332449, 944.028, -25061.5, -0.000299453, 0.00332449, -0.000299453, -13404.8])
  #stress_alpha = np.array([-18284.7, -1808.06, 0.00167735, -1808.06, -12392.2, 0.00069253, 0.00167735, 0.00069253, -13989.4])
  #p = calcData(p, df_dsigma, P_n, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.0723199, -0.906983, 0.326434, 0.256111, 1.53284e-07, -9.23145e-08])
  #P_n = np.array([-0.274192, -0.938804, 0.0433213, 0.203932, 1.22055e-07, -7.35067e-08])
  #alpha = 42353.5
  #s0     = np.array([-14486.9, 48.3381, 0.00199985, 48.3381, -14644.4, 1.1731e-05, 0.00199985, 1.1731e-05, -9333.14])
  #strial = np.array([-18235.9, 496.183, 0.00266217, 496.183, -19853, -0.000143861, 0.00266217, -0.000143861, -11369])
  #salpha = np.array([-6622.89, -8141.06, -0.00250727, -8141.06, 19908.7, 0.00296941, -0.00250727, 0.00296941, -13203.8])
  #p = calcData(p, df_dsigma, P_n, s0, strial, salpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.0723199, -0.906983, 0.326434, 0.256111, 1.53284e-07, -9.23145e-08])
  #P_n = np.array([-0.274192, -0.938804, 0.0433213, 0.203932, 1.22055e-07, -7.35067e-08])
  #s0     = np.array([-14486.9, 48.3381, 0.00199985, 48.3381, -14644.4, 1.1731e-05, 0.00199985, 1.1731e-05, -9333.14])
  #strial = np.array([-16361.4, 272.261, 0.00233101, 272.261, -17248.7, -6.60649e-05, 0.00233101, -6.60649e-05, -10351])
  #salpha = np.array([-15436.3, -415.762, 0.00191923, -415.762, -14081.3, 0.000181931, 0.00191923, 0.000181931, -10497.2])
  #p = calcData(p, df_dsigma, P_n, s0, strial, salpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.906983, -0.0723199, 0.326434, -0.256111, 4.80152e-07, 1.01684e-07])
  #Pn = np.array([-0.938804, -0.274192, 0.0433213, -0.203932, 3.82328e-07, 8.09677e-08])
  #stress_old   = np.array([-15229.8, -192.864, 0.00156048, -192.864, -14601.2, 9.43078e-05, 0.00156048, 9.43078e-05, -11211.8])
  #stress_trial = np.array([-18978.8, 254.981, 0.0022228, 254.981, -19809.8, -6.12842e-05, 0.0022228, -6.12842e-05, -13247.6])
  #stress_alpha = np.array([-12644.2, 1631.03, -0.000356982, 1631.03, -17959.7, -0.000607619, -0.000356982, -0.000607619, -13539.9])
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.0723199058966196, -0.906982710164428, 0.326434205353643, 0.256111116996115, 1.53283933944362e-07, -9.23144656544059e-08])
  #stress_old = np.array([-14433.5, 295.662, 0.00102008, 295.662, -15397.1, -9.40987e-05, 0.00102008, -9.40987e-05, -11800.7])
  #stress_trial = np.array([-15370.8, 407.624, 0.00118567, 407.624, -16699.2, -0.000132997, 0.00118567, -0.000132997, -12309.7])
  #stress_alpha = np.array([-14908.3, 63.6121, 0.000979773, 63.6121, -15115.6, -8.99886e-06, 0.000979773, -8.99886e-06, -12382.7])
  #Pn = np.array([-0.274192,    -0.938804,    0.0433213,     0.203932,  1.22054e-07, -7.35067e-08])
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.0723199058938287, -0.906982710166111, 0.326434205353274, 0.256111116991412, 1.5328393384551e-07, -9.23144663005642e-08])
  #stress_old = np.array([-15055, 172.774, 0.00104511, 172.774, -15618.1, -4.83461e-05, 0.00104511, -4.83461e-05, -12359.5])
  #stress_trial = np.array([-16929.5, 396.697, 0.00137627, 396.697, -18222.4, -0.000126142, 0.00137627, -0.000126142, -13377.5])
  #stress_alpha = np.array([-16004.5, -291.326, 0.000964483, -291.326, -15055.1, 0.000121854, 0.000964483, 0.000121854, -13523.6])
  #Pn = np.array([-0.274192, -0.938804, 0.0433213, 0.203932, 1.22054e-07, -7.35067e-08])
  #Pn_calc = computePn(df_dsigma)
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)
  #print("Pn = ", Pn, "Pn_calc = ", Pn_calc)

  #df_dsigma = np.array([-0.906983, -0.0723199, 0.326434, -0.256111, 4.80152e-07, 1.01684e-07])
  #P_n = np.array([-0.938804, -0.274192, 0.0433213, -0.203932, 3.82328e-07, 8.09677e-08])
  #s0 = np.array([-15966.6, -319.478, 0.000947633, -319.478, -14925.5, 0.000132001, 0.000947633, 0.000132001, -13529.6])
  #strial = np.array([-23464.7, 576.211, 0.00227228, 576.211, -25342.6, -0.000179183, 0.00227228, -0.000179183, -17601.2])
  #salpha = np.array([56058.6, 17850.7, -0.0301136, 17850.7, -2116.58, -0.00703771, -0.0301136, -0.00703771, -21270.9])
  #p = calcData(p, df_dsigma, Pn, s0, strial, salpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.906982710165378, -0.0723199058949562, 0.326434205353318, -0.256111116993226, 4.801524242618e-07, 1.01684448941657e-07])
  #stress_old = np.array([-15966.6, -319.478, 0.000947633, -319.478, -14925.5, 0.000132001, 0.000947633, 0.000132001, -13529.6])
  #stress_trial = np.array([-23464.7, 576.211, 0.00227228, 576.211, -25342.6, -0.000179183, 0.00227228, -0.000179183, -17601.2])
  #stress_alpha = np.array([-10795.4, 3328.3, -0.00288729, 3328.3, -21642.3, -0.00127185, -0.00288729, -0.00127185, -18185.9])
  #Pn = np.array([-0.938804, -0.274192, 0.0433213, -0.203932, 3.82328e-07, 8.09677e-08])
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.906983, -0.0723199, 0.326434, -0.256111, 4.80152e-07, 1.01684e-07])
  #Pn = np.array([-0.938804, -0.274192, 0.0433213, -0.203932, 3.82328e-07, 8.09677e-08])
  #s0 = np.array([-15966.6, -319.478, 0.000947633, -319.478, -14925.5, 0.000132001, 0.000947633, 0.000132001, -13529.6])
  #strial = np.array([-19715.7, 128.366, 0.00160996, 128.366, -20134, -2.35908e-05, 0.00160996, -2.35908e-05, -15565.4])
  #salpha = np.array([20046, 8765.61, -0.014583, 8765.61, -8521.01, -0.00345286, -0.014583, -0.00345286, -17400.2])
  #p = calcData(p, df_dsigma, Pn, s0, strial, salpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.906982710165378, -0.0723199058949562, 0.326434205353318, -0.256111116993226, 4.801524242618e-07, 1.01684448941657e-07])
  #stress_old = np.array([-15966.6, -319.478, 0.000947633, -319.478, -14925.5, 0.000132001, 0.000947633, 0.000132001, -13529.6])
  #stress_trial = np.array([-19715.7, 128.366, 0.00160996, 128.366, -20134, -2.35908e-05, 0.00160996, -2.35908e-05, -15565.4])
  #stress_alpha = np.array([-13381, 1504.41, -0.00096983, 1504.41, -18283.9, -0.000569926, -0.00096983, -0.000569926, -15857.7])
  #Pn = np.array([-0.938804, -0.274192, 0.0433213, -0.203932, 3.82328e-07, 8.09677e-08])
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.906983, -0.0723199, 0.326434, -0.256111, 4.80152e-07, 1.01684e-07])
  #Pn = np.array([-0.938804, -0.274192, 0.0433213, -0.203932, 3.82328e-07, 8.09677e-08])
  #s0 = np.array([-15966.6, -319.478, 0.000947633, -319.478, -14925.5, 0.000132001, 0.000947633, 0.000132001, -13529.6])
  #strial = np.array([-17841.2, -95.556, 0.00127879, -95.556, -17529.7, 5.42052e-05, 0.00127879, 5.42052e-05, -14547.5])
  #salpha = np.array([2039.68, 4223.06, -0.00681768, 4223.06, -11723.2, -0.00166043, -0.00681768, -0.00166043, -15464.9])
  #p = calcData(p, df_dsigma, Pn, s0, strial, salpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.906982710165378, -0.0723199058949562, 0.326434205353318, -0.256111116993226, 4.801524242618e-07, 1.01684448941657e-07])
  #stress_old = np.array([-15966.6, -319.478, 0.000947633, -319.478, -14925.5, 0.000132001, 0.000947633, 0.000132001, -13529.6])
  #stress_trial = np.array([-17841.2, -95.556, 0.00127879, -95.556, -17529.7, 5.42052e-05, 0.00127879, 5.42052e-05, -14547.5])
  #stress_alpha = np.array([-14673.8, 592.467, -1.10982e-05, 592.467, -16604.7, -0.000218962, -1.10982e-05, -0.000218962, -14693.7])
  #Pn = np.array([-0.938804, -0.274192, 0.0433213, -0.203932, 3.82328e-07, 8.09677e-08])
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.906983, -0.0723199, 0.326434, -0.256111, 4.80152e-07, 1.01684e-07])
  #Pn = np.array([-0.938804, -0.274192, 0.0433213, -0.203932, 3.82328e-07, 8.09677e-08])
  #s0 = np.array([-15966.6, -319.478, 0.000947633, -319.478, -14925.5, 0.000132001, 0.000947633, 0.000132001, -13529.6])
  #strial = np.array([-16903.9, -207.517, 0.00111321, -207.517, -16227.6, 9.31031e-05, 0.00111321, 9.31031e-05, -14038.6])
  #salpha = np.array([-6963.48, 1951.79, -0.00293502, 1951.79, -13324.3, -0.000764213, -0.00293502, -0.000764213, -14497.3])
  #p = calcData(p, df_dsigma, Pn, s0, strial, salpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.906982710165378, -0.0723199058949562, 0.326434205353318, -0.256111116993226, 4.801524242618e-07, 1.01684448941657e-07])
  #stress_old = np.array([-15966.6, -319.478, 0.000947633, -319.478, -14925.5, 0.000132001, 0.000947633, 0.000132001, -13529.6])
  #stress_trial = np.array([-16903.9, -207.517, 0.00111321, -207.517, -16227.6, 9.31031e-05, 0.00111321, 9.31031e-05, -14038.6])
  #Pn = np.array([-0.938804, -0.274192, 0.0433213, -0.203932, 3.82328e-07, 8.09677e-08])
  #stress_alpha = np.array([-16220.3679148487, -59.0384150464422, 0.000834848486827511, -59.0384150464422, -16027.9624869284, 3.41522281617778e-05, 0.000834848486827511, 3.41522281617778e-05, -14070.0929779496])
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)

  #dg_dsigma = np.array([-0.0723199527894327, -0.90698268187496, 0.326434211554792, 0.25611119603414, 3.65961078888307e-08, 4.7647574260901e-08])
  #stress_old = np.array([-15426.9, 180.083, 0.000240119, 180.083, -16013.8, 0.000139581, 0.000240119, 0.000139581, -13552.9])
  #stress_trial = np.array([-16364.2, 292.045, 0.000290286, 292.045, -17315.9, 0.000177317, 0.000290286, 0.000177317, -14061.9])
  #Pn = np.array([-0.274192, -0.938804, 0.0433213, 0.203932, 2.91402e-08, 3.79401e-08])
  #stress_alpha = np.array([-16048.3904500766, 57.195515415857, 0.000256728448292808, 57.195515415857, -16234.7898176814, 0.000133625503330731, 0.000256728448292808, 0.000133625503330731, -14111.7424688338])
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.906982681855192, -0.072319952822187, 0.326434211559107, -0.256111196089332, 1.82978512696659e-07, 6.09973706483911e-08])
  #stress_old = np.array([-16371.7, -287.437, 0.000241652, -287.437, -15435, 8.64158e-05, 0.000241652, 8.64158e-05, -14717.7])
  #stress_trial = np.array([-17309, -175.476, 0.000291819, -175.476, -16737.1, 0.000124152, 0.000291819, 0.000124152, -15226.6])
  #Pn = np.array([-0.938804, -0.274192, 0.0433213, -0.203932, 1.45699e-07, 4.857e-08])
  #stress_alpha = np.array([-16625.4385332545, -26.9862662168994, 0.000185730972597737, -26.9862662168994, -16537.4906772603, 8.87865962175641e-05, 0.000185730972597737, 8.87865962175641e-05, -15258.1912724386])
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.906984031968465, -0.0723177148639539, 0.326433915610787, -0.256107423965131, 1.82978400969841e-07, 6.09966938022202e-08])
  #stress_old = np.array([-17312.7, -2.97641, 1.9429e-06, -2.97641, -17303, 6.18031e-07, 1.9429e-06, 6.18031e-07, -17299.8])
  #stress_trial = np.array([-17320.1, -2.10171, 2.33483e-06, -2.10171, -17313.2, 9.12844e-07, 2.33483e-06, 9.12844e-07, -17303.8])
  #Pn = np.array([-0.938805, -0.27419, 0.0433213, -0.203929, 1.45699e-07, 4.85695e-08])
  #stress_alpha = np.array([-17314.7557326165, -0.946900118809321, 1.50976740492499e-06, -0.946900118809321, -17311.6696856996, 6.37805781679766e-07, 1.50976740492499e-06, 6.37805781679766e-07, -17304.0016537346])
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)

  #df_dsigma = np.array([-0.114273679868487, -0.869617765587704, 0.310261132066981, 0.336639775685172, -0.108970099213085, 0.0961440413973013])
  #stress_old = np.array([-17187.2, 21.3857, -19.9296, 21.3857, -17238.6, 7.03823, -19.9296, 7.03823, -17106.7])
  #stress_trial = np.array([-17313.6, 39.901, -31.2799, 39.901, -17409.6, 12.735, -31.2799, 12.735, -17187.7])
  #Pn = np.array([-0.310833, -0.905299, 0.0232813, 0.26494, -0.0857609, 0.0756667])
  #stress_alpha = np.array([-16325.3518481194, 172.098900885759, -122.662126341232, 172.098900885759, -16744.4391009337, 54.1794689965474, -122.662126341232, 54.1794689965474, -15832.6812521302])
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)
  #Pn_calc = computePn(df_dsigma)
  #print(np.linalg.norm(stress_old - stress_trial))
  #print("I1, J2 = ", I1J2(stress_alpha.reshape(3,3)))
  #print("Pn_calc = ", Pn_calc)
 

  #df_dsigma = np.array([-0.107495852042674, -0.874137273061191, 0.309817061590607, 0.328594761931334, -0.107077794153276, 0.0943498087792677])
  #stress_old = np.array([-15042, 405.896, -262.725, 405.896, -16023.4, 125.881, -262.725, 125.881, -13989.4])
  #stress_trial = np.array([-26831.1, 2133.6, -1321.85, 2133.6, -31974.3, 657.461, -1321.85, 657.461, -21541.8])
  #stress_alpha = np.array([-20321.7, -3387.4, 477.253, -3387.4, -12583.9, -927.79, 477.253, -927.79, -22044.1])
  #Pn = np.array([-0.305217, -0.909193, 0.0235509, 0.258874, -0.0843581, 0.0743307])
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)
  #Pn_calc = computePn(df_dsigma)
  #print(np.linalg.norm(stress_old - stress_trial))
  #print("I1, J2 = ", I1J2(stress_alpha.reshape(3,3)))
  #print("Pn_calc = ", Pn_calc)

  #df_dsigma = np.array([-0.577350269189626, -0.577350269189626, -0.577350269189626, 0, 0, 0])
  #stress_old = np.array([ -17980.9, 96.9065, -60.5767, 96.9065, -18214.7, 29.9, -60.5767, 29.9, -17738.5])
  #stress_trial = np.array([-19307.2, 291.273, -179.728, 291.273, -20009.1, 89.7027, -179.728, 89.7027, -18588.1])
  #Pn = np.array([-0.57735, -0.57735, -0.57735, 0, 0, 0])
  #stress_alpha = np.array([-14002.4599769023, 595.271579525215, -376.176325618175, 595.271579525215, -15439.3486369625, 183.959027079059, -376.176325618175, 183.959027079059, -12496.3825238968])
  #p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)
  #Pn_calc = computePn(df_dsigma)
  #print("I1, J2 = ", I1J2(stress_trial.reshape(3,3)))
  #print("Pn_calc = ", Pn_calc)

  p.show_grid()
  p.show()

  return

plotData()
