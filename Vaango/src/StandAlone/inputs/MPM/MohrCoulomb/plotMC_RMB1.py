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
  #if (J3 < 0) :
  #  q = -q
  r_cubed = 27.0/2.0 * J3
  cos3theta = r_cubed / (q * q * q)
  if (cos3theta < -1):
    cos3theta = -1
  if (cos3theta > 1):
    cos3theta = 1
  theta = np.arccos(cos3theta) / 3.0

  xi = 1.0 / np.sqrt(3.0) * I1
  rho = np.sqrt(2.0 / 3.0) * q

  sigma_1 = p + 2.0 / 3.0 * q * np.cos(theta)
  sigma_2 = p + 2.0 / 3.0 * q * np.cos(theta - 2.0 * np.pi / 3.0)
  sigma_3 = p + 2.0 / 3.0 * q * np.cos(theta + 2.0 * np.pi / 3.0)

  sig_prin = np.array([sigma_1, sigma_2, sigma_3])

  return sig_prin

def mohr_coulomb(c, phi, low, high, numcells):

    rad = np.linspace(0.01, high, numcells)
    grid = pv.CylinderStructured(radius = rad, height = np.sqrt(3)*(high - low),
                                 center = [0.5*(high+low), 0.5*(high+low), 0.5*(high+low)],
                                 direction = [1, 1, 1],
                                 theta_resolution = numcells,
                                 z_resolution = np.int(1.73*numcells))

    vtkGrid = pv.UnstructuredGrid(grid)
    stresses = vtkGrid.points

    sin_phi = np.sin(phi * np.pi / 180.0)
    cos_phi = np.cos(phi * np.pi / 180.0)

    yield_fn = []
    for stress in stresses:
      ind = np.argsort(stress)
      i1 = ind[2]
      i3 = ind[0]
      f = ((stress[i1] - stress[i3]) -
           (stress[i1] + stress[i3]) * sin_phi - 
           2.0 * c * cos_phi) / \
           (abs(stress[i1]) + abs(stress[i3]) + 2.0 * c)
      yield_fn.append(f)

    vtkGrid["data"] = np.array(yield_fn)
    return vtkGrid

def drawArrow(start, end, scale):

  length = np.linalg.norm(end - start)
  direction = (end - start) / length
  cyl_bot = start
  cyl_top = end
  cyl_cen = 0.5 * (cyl_bot + cyl_top)
  arrow_cyl = pv.Cylinder(cyl_cen, direction = direction, radius = scale, 
                          height = length)
  cone_bot = cyl_top
  cone_top = cone_bot + direction * length * 0.3
  cone_cen = 0.5 * (cone_bot + cone_top)
  arrow_cone = pv.Cone(cone_cen, direction, height = np.linalg.norm(cone_top - cone_bot), 
                       radius = 2.5*scale)

  return arrow_cyl, arrow_cone
  
def toMatrix(vec):
  return np.array([vec[0], vec[3], vec[4], vec[3], vec[1], vec[5], vec[4], vec[5], vec[2]]).reshape(3,3)

def rotateToInit(eigvec, stress):
  return np.matmul(np.matmul(eigvec.transpose(), stress), eigvec)

def calcData(p, df_dsigma, P_n, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad):

  df_dsigma = toMatrix(df_dsigma)
  P_n = toMatrix(P_n)

  alpha = np.linalg.norm(stress_trial - stress_alpha)

  stress_old = stress_old.reshape(3,3)
  stress_trial = stress_trial.reshape(3,3)
  stress_alpha = stress_alpha.reshape(3,3)

  stress_old_eig = principal_stresses(stress_old)
  stress_trial_eig = principal_stresses(stress_trial)
  stress_alpha_eig = principal_stresses(stress_alpha)

  print(stress_old_eig)
  print(stress_trial_eig)
  print(stress_alpha_eig)

  stress_old_ball = pv.Sphere(radius=ball_rad, center = stress_old_eig)
  stress_trial_ball = pv.Sphere(radius=ball_rad, center = stress_trial_eig)
  stress_alpha_ball = pv.Sphere(radius=ball_rad, center = stress_alpha_eig)

  stress_Pn = stress_trial - P_n * alpha 
  stress_df = stress_old + df_dsigma * arrow_len
  stress_Pn_eig = principal_stresses(stress_Pn)
  stress_df_eig = principal_stresses(stress_df)

  stress_Pn_ball = pv.Sphere(radius=ball_rad, center = stress_Pn_eig)
  stress_df_ball = pv.Sphere(radius=ball_rad, center = stress_df_eig)

  df_dsigma_cyl, df_dsigma_cone = drawArrow(stress_old_eig, stress_df_eig, cyl_rad)
  P_n_cyl, P_n_cone = drawArrow(stress_trial_eig, stress_Pn_eig, cyl_rad)
  sigma_cyl, sigma_cone = drawArrow(stress_old_eig, stress_trial_eig, cyl_rad)

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
  grid = mohr_coulomb(c, phi, -1.0e5, upper_lim, 100)
  contours = grid.contour([0.0])

  pv.set_plot_theme('document')
  p = pv.Plotter(notebook=False)
  #p.add_mesh(grid, opacity=0.3)
  p.add_mesh(contours, color="#FAEBD7", opacity=1.0, show_scalar_bar=False)
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

  df_dsigma = np.array([-0.0723199, -0.906983, 0.326434, 0.256111, 1.53284e-07, -9.23145e-08])
  P_n = np.array([-0.274192, -0.938804, 0.0433213, 0.203932, 1.22055e-07, -7.35067e-08])
  s0     = np.array([-14486.9, 48.3381, 0.00199985, 48.3381, -14644.4, 1.1731e-05, 0.00199985, 1.1731e-05, -9333.14])
  strial = np.array([-16361.4, 272.261, 0.00233101, 272.261, -17248.7, -6.60649e-05, 0.00233101, -6.60649e-05, -10351])
  salpha = np.array([-15436.3, -415.762, 0.00191923, -415.762, -14081.3, 0.000181931, 0.00191923, 0.000181931, -10497.2])
  p = calcData(p, df_dsigma, P_n, s0, strial, salpha, ball_rad, arrow_len, cyl_rad)

  df_dsigma = np.array([-0.906983, -0.0723199, 0.326434, -0.256111, 4.80152e-07, 1.01684e-07])
  Pn = np.array([-0.938804, -0.274192, 0.0433213, -0.203932, 3.82328e-07, 8.09677e-08])
  stress_old   = np.array([-15229.8, -192.864, 0.00156048, -192.864, -14601.2, 9.43078e-05, 0.00156048, 9.43078e-05, -11211.8])
  stress_trial = np.array([-18978.8, 254.981, 0.0022228, 254.981, -19809.8, -6.12842e-05, 0.0022228, -6.12842e-05, -13247.6])
  stress_alpha = np.array([-12644.2, 1631.03, -0.000356982, 1631.03, -17959.7, -0.000607619, -0.000356982, -0.000607619, -13539.9])
  p = calcData(p, df_dsigma, Pn, stress_old, stress_trial, stress_alpha, ball_rad, arrow_len, cyl_rad)

  p.show_grid()
  p.show()

  return

plotData()
