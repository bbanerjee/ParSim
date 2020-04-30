import numpy as np
import pyvista as pv

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

def J2(s):

  J2 = ((s[0,0] - s[1,1]) * (s[0,0] - s[1,1]) + (s[0,0] - s[2,2]) * (s[0,0] - s[2,2]) + \
                 (s[1,1] - s[2,2]) * (s[1,1] - s[2,2])) / 6.0 + \
       (s[0,1] * s[0,1] + s[0,2] * s[0,2] + s[1,2] * s[1,2])
  return J2

def J3(s):
  I1 = s[0,0] + s[1,1] + s[2,2]
  I2 = s[0,0] * s[1,1] + s[1,1] * s[2,2] + s[2,2] * s[0,0] - s[0,1] * s[0,1] - s[0,2] * s[0,2] - s[1,2] * s[1,2]
  I3 = s[0][0] * s[1,1] * s[2,2] + 2 * s[0,1] * s[0,2] * s[1,2] - s[0,0] * s[1,2] * s[1,2] - s[1,1] * s[0,2] * s[0,2] - s[2,2] * s[0,1] * s[0,1]
  J3 = I1 * I1 * I1 * 2.0 / 27.0 - I1 * I2 / 3.0 + I3
  return J3

def lodeAngle(s):

  cos_theta = J3(s)/2 * (3/J2(s))**(1.5)
  return np.arccos(cos_theta) * 180.0 /np.pi

def drawArrow(start, direction, length, scale):

  cyl_bot = start
  cyl_top = start + direction * length
  cyl_cen = 0.5 * (cyl_bot + cyl_top)
  arrow_cyl = pv.Cylinder(cyl_cen, direction = direction, radius = scale, 
                          height = length)
  cone_bot = cyl_top
  cone_top = cone_bot + direction * length * 0.3
  cone_cen = 0.5 * (cone_bot + cone_top)
  arrow_cone = pv.Cone(cone_cen, direction, height = np.linalg.norm(cone_top - cone_bot), 
                       radius = 2.5*scale)

  return arrow_cyl, arrow_cone
  
def plotData1():
  df_dsigma = np.array([0.316228, 0, -0.948683])

  P_n = np.array([6.32456e+06, -3.16228e+07, -1.45465e+08])
  P_n_length = np.linalg.norm(P_n);
  P_n = P_n / P_n_length;

  stress_init = \
    np.array([294450, 1.28206e-12, 1.2255e-12, 
              1.28206e-12, 86602.9, -5.02731e-32, 
              1.2255e-12, -5.02731e-32, 86602.9])
  stress_trial = \
    np.array([299422, 1.31845e-12, 1.31741e-12, 
              1.31845e-12, 88065.2, -5.16843e-32, 
              1.31741e-12, -5.16843e-32, 88065.2])
  stress_closest = \
    np.array([298836, 1.31845e-12, 1.31741e-12, 
              1.31845e-12, 88065.2, -5.16843e-32, 
              1.31741e-12, -5.16843e-32, 89821.3])
  stress_new = \
    np.array([299338, 1.31845e-12, 1.31741e-12, 
              1.31845e-12, 88483.3, -5.16843e-32, 
              1.31741e-12, -5.16843e-32, 89988.5])

  stress_init_eig, vec = np.linalg.eigh(np.reshape(stress_init, (3,3)))
  stress_trial_eig, vec = np.linalg.eigh(np.reshape(stress_trial, (3,3)))
  stress_closest_eig, vec = np.linalg.eigh(np.reshape(stress_closest, (3,3)))
  stress_new_eig, vec = np.linalg.eigh(np.reshape(stress_new, (3,3)))
  stress_init_eig = np.flip(stress_init_eig)
  stress_trial_eig = np.flip(stress_trial_eig)
  stress_closest_eig = np.flip(stress_closest_eig)
  stress_new_eig = np.flip(stress_new_eig)

  stress_init_ball = pv.Sphere(radius=1000, center = stress_init_eig)
  stress_trial_ball = pv.Sphere(radius=1000, center = stress_trial_eig)
  stress_closest_ball = pv.Sphere(radius=1000, center = stress_closest_eig)
  stress_new_ball = pv.Sphere(radius=1000, center = stress_new_eig)

  df_dsigma_cyl, df_dsigma_cone = drawArrow(stress_init_eig, df_dsigma, 1.0e5, 1000)
  P_n_cyl, P_n_cone = drawArrow(stress_trial_eig, P_n, 1.0e5, 1000)

  stress_dir = stress_trial_eig - stress_init_eig
  sigma_cyl, sigma_cone = drawArrow(stress_init_eig, 
                                    stress_dir/np.linalg.norm(stress_dir), 
                                    1.0e5, 1000)

  c = 1.0e4
  phi = 30
  grid = mohr_coulomb(c, phi, -1.0e5, 5.0e5, 100)
  contours = grid.contour([0.0])

  blue = np.array([12/256, 238/256, 246/256])
  grey = np.array([189/256, 189/256, 189/256])
  yellow = np.array([255/256, 247/256, 0/256])

  pv.set_plot_theme('document')
  p = pv.Plotter(notebook=False)
  #p.add_mesh(grid, opacity=0.3)
  p.add_mesh(contours, opacity=0.6, show_scalar_bar=False)
  p.add_mesh(stress_init_ball, color="red", show_edges=False)
  p.add_mesh(stress_trial_ball, color="red", show_edges=False)
  p.add_mesh(stress_closest_ball, color=blue, show_edges=False)
  p.add_mesh(stress_new_ball, color=yellow, show_edges=False)
  p.add_mesh(df_dsigma_cyl, color="green", show_edges=False)
  p.add_mesh(df_dsigma_cone, color="green", show_edges=False)
  p.add_mesh(P_n_cyl, color=yellow, show_edges=False)
  p.add_mesh(P_n_cone, color=yellow, show_edges=False)
  p.add_mesh(sigma_cyl, color="blue", show_edges=False)
  p.add_mesh(sigma_cone, color="blue", show_edges=False)
  #p.add_axes()
  #p.add_bounding_box()
  p.show_grid()
  p.show()

  return

def plotData2():

  fac = 0.1
  upper_lim = 1.0e5 * fac
  ball_rad = 2000 * fac
  cyl_rad = 1000 * fac
  arrow_len = 1.0e5 * fac 

  df_dsigma = np.array([0.316228, 0, -0.948683])

  P_n = np.array([0.0424476, -0.212238, -0.976296])
  P_n_length = np.linalg.norm(P_n);
  P_n = P_n / P_n_length;

  alpha = 30226

  stress_init = \
    np.array([-5043.43, 0, 0, 0, -8879.05, 0, 0, 0, -13228.2]).reshape(3,3)
  stress_trial = \
    np.array([-13336.4, -1.59465e-10, 6.57669e-10, -1.59465e-10, -23479, -1.92522e-09, 6.57669e-10, -1.92522e-09, -34979.4]).reshape(3,3)
  stress_alpha = \
    np.array([-14502.7, -1.59419e-10, 6.57582e-10, -1.59419e-10, -17647, -1.92417e-09, 6.57582e-10, -1.92417e-09, -8152.65]).reshape(3,3)

  stress_init_eig = np.diagonal(stress_init)
  stress_trial_eig = np.diagonal(stress_trial)
  stress_alpha_eig = np.diagonal(stress_alpha)

  # The eigenvectors of P_n are not the same as those of stress
  stress_alpha_ball = pv.Sphere(radius=ball_rad, center = stress_alpha_eig)
  stress_alpha_actual = stress_trial_eig - P_n * alpha

  #print(stress_alpha_eig)
  #print(stress_alpha_actual)

  stress_init_ball = pv.Sphere(radius=ball_rad, center = stress_init_eig)
  stress_trial_ball = pv.Sphere(radius=ball_rad, center = stress_trial_eig)
  stress_actual_ball = pv.Sphere(radius=ball_rad, center = stress_alpha_actual)

  df_dsigma_cyl, df_dsigma_cone = drawArrow(stress_init_eig, df_dsigma, arrow_len, cyl_rad)
  P_n_cyl, P_n_cone = drawArrow(stress_trial_eig, -P_n, arrow_len, cyl_rad)

  stress_dir = stress_trial_eig - stress_init_eig
  sigma_cyl, sigma_cone = drawArrow(stress_init_eig, 
                                    stress_dir/np.linalg.norm(stress_dir), 
                                    arrow_len, cyl_rad)

  strial = np.array([-9189.88, -7.97093e-11, 3.28791e-10, -7.97093e-11, -16179, -9.62086e-10, 3.28791e-10, -9.62086e-10, -24103.6]).reshape(3,3)
  snew   = np.array([-9596.73, -7.97093e-11, 3.28791e-10, -7.97093e-11, -14144.7, -9.62086e-10, 3.28791e-10, -9.62086e-10, -14745.9]).reshape(3,3)

  strial_eig = np.diagonal(strial)
  snew_eig = np.diagonal(snew)

  strial_ball = pv.Sphere(radius=ball_rad, center = strial_eig)
  snew_ball = pv.Sphere(radius=ball_rad, center = snew_eig)

  s_dir = strial_eig - stress_init_eig
  s_cyl, s_cone = drawArrow(stress_init_eig, 
                                    s_dir/np.linalg.norm(s_dir), 
                                    arrow_len, cyl_rad)
  sP_n_cyl, sP_n_cone = drawArrow(strial_eig, -P_n, arrow_len, cyl_rad)

  strial2 = np.array([-13743.2, -1.59419e-10, 6.57582e-10, -1.59419e-10, -21444.6, -1.92417e-09, 6.57582e-10, -1.92417e-09, -25621.4]).reshape(3,3)
  strial2_eig = np.diagonal(strial2)
  strial2_ball = pv.Sphere(radius=ball_rad, center = strial2_eig)
  s2_dir = strial2_eig - snew_eig
  s2_cyl, s2_cone = drawArrow(snew_eig, s2_dir/np.linalg.norm(s2_dir), arrow_len, cyl_rad)
  s2P_n_cyl, s2P_n_cone = drawArrow(strial2_eig, -P_n, arrow_len, cyl_rad)

  strial3 = np.array([-11670, -1.19564e-10, 4.93187e-10, -1.19564e-10, -17794.6, -1.44313e-09, 4.93187e-10, -1.44313e-09, -20183.7]).reshape(3,3)
  strial3_eig = np.diagonal(strial3)
  strial3_ball = pv.Sphere(radius=ball_rad, center = strial3_eig)
  s3_dir = strial3_eig - snew_eig
  s3_cyl, s3_cone = drawArrow(snew_eig, s3_dir/np.linalg.norm(s3_dir), arrow_len, cyl_rad)
  s3P_n_cyl, s3P_n_cone = drawArrow(strial3_eig, -P_n, arrow_len, cyl_rad)

  strial4 = np.array([-10633.3, -9.96366e-11, 4.10989e-10, -9.96366e-11, -15969.7, -1.20261e-09, 4.10989e-10, -1.20261e-09, -17464.8]).reshape(3,3)
  snew4 = np.array([-10798, -9.96366e-11, 4.10989e-10, -9.96366e-11, -15146.3, -1.20261e-09, 4.10989e-10, -1.20261e-09, -13677.5]).reshape(3,3)
  strial4_eig = np.diagonal(strial4)
  snew4_eig = np.diagonal(snew4)
  strial4_ball = pv.Sphere(radius=ball_rad, center = strial4_eig)
  snew4_ball = pv.Sphere(radius=ball_rad, center = snew4_eig)
  s4_dir = strial4_eig - snew_eig
  s4_cyl, s4_cone = drawArrow(snew_eig, s4_dir/np.linalg.norm(s4_dir), arrow_len, cyl_rad)
  s4P_n_cyl, s4P_n_cone = drawArrow(strial4_eig, -P_n, arrow_len, cyl_rad)

  c = 1.0e4
  phi = 30
  grid = mohr_coulomb(c, phi, -1.0e5, upper_lim, 100)
  contours = grid.contour([0.0])

  blue = np.array([12/256, 238/256, 246/256])
  grey = np.array([189/256, 189/256, 189/256])
  yellow = np.array([255/256, 247/256, 0/256])

  pv.set_plot_theme('document')
  p = pv.Plotter(notebook=False)
  #p.add_mesh(grid, opacity=0.3)
  p.add_mesh(contours, opacity=0.6, show_scalar_bar=False)
  #p.add_mesh(stress_init_ball, color="red", show_edges=False)
  #p.add_mesh(stress_trial_ball, color=blue, show_edges=False)
  #p.add_mesh(stress_alpha_ball, color="#F4C2C2", show_edges=False)
  #p.add_mesh(stress_actual_ball, color=yellow, show_edges=False)
  p.add_mesh(df_dsigma_cyl, color="green", show_edges=False)
  p.add_mesh(df_dsigma_cone, color="green", show_edges=False)
  #p.add_mesh(P_n_cyl, color=yellow, show_edges=False)
  #p.add_mesh(P_n_cone, color=yellow, show_edges=False)
  #p.add_mesh(sigma_cyl, color="blue", show_edges=False)
  #p.add_mesh(sigma_cone, color="blue", show_edges=False)

  #p.add_mesh(strial_ball, color=grey, show_edges=False)
  #p.add_mesh(snew_ball, color="#ff7E00", show_edges=False)
  #p.add_mesh(s_cyl, color="#ff7E00", show_edges=False)
  #p.add_mesh(s_cone, color="#ff7E00", show_edges=False)
  #p.add_mesh(sP_n_cyl, color=yellow, show_edges=False)
  #p.add_mesh(sP_n_cone, color=yellow, show_edges=False)

  #p.add_mesh(strial2_ball, color=grey, show_edges=False)
  #p.add_mesh(s2_cyl, color="#ff7E00", show_edges=False)
  #p.add_mesh(s2_cone, color="#ff7E00", show_edges=False)
  #p.add_mesh(s2P_n_cyl, color=yellow, show_edges=False)
  #p.add_mesh(s2P_n_cone, color=yellow, show_edges=False)

  #p.add_mesh(strial3_ball, color=grey, show_edges=False)
  #p.add_mesh(s3_cyl, color="#ff7E00", show_edges=False)
  #p.add_mesh(s3_cone, color="#ff7E00", show_edges=False)
  #p.add_mesh(s3P_n_cyl, color=yellow, show_edges=False)
  #p.add_mesh(s3P_n_cone, color=yellow, show_edges=False)

  p.add_mesh(strial4_ball, color=grey, show_edges=False)
  p.add_mesh(snew4_ball, color="#ff7E00", show_edges=False)
  p.add_mesh(s4_cyl, color="#ff7E00", show_edges=False)
  p.add_mesh(s4_cone, color="#ff7E00", show_edges=False)
  p.add_mesh(s4P_n_cyl, color=yellow, show_edges=False)
  p.add_mesh(s4P_n_cone, color=yellow, show_edges=False)

  print(J2(strial4), J2(snew4))
  print(J3(strial4), J3(snew4))
  print(lodeAngle(strial4), lodeAngle(snew4))

  #p.add_axes()
  #p.add_bounding_box()
  p.show_grid()
  p.show()

  return

plotData2()
