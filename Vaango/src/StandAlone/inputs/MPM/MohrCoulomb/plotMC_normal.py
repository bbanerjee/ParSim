import numpy as np
import pyvista as pv

def I1J2J3(stress):

  I = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).reshape(3,3)
  I1 = stress.trace()
  s = stress - I * (I1 / 3.0)
  ss = np.matmul(s, s)
  J2 = ss.trace() / 2.0
  J3 = np.linalg.det(s)
  return I1, J2, J3
  
def pqtheta(stress):

  I1, J2, J3 = I1J2J3(stress)

  p = I1 / 3.0
  q = np.sqrt(3.0 * J2)
  r_cubed = 27.0/2.0 * J3
  cos3theta_c = r_cubed / (q * q * q)
  if (cos3theta_c < -1):
    cos3theta_c = -1
  if (cos3theta_c > 1):
    cos3theta_c = 1
  theta_c = np.arccos(cos3theta_c) / 3.0

  return p, q, theta_c

def df_dsigma_fn(stress, phi):
  sin_phi = np.sin(phi * np.pi/180)
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

  R, dR_dtheta = dR_dtheta_fn(phi, theta)
  dq_dsigma, dtheta_dsigma = dtheta_dsigma_fn(s, ss, J2, q, J3, theta)
  dp_dsigma = dp_dsigma_fn()
  dR_dsigma = dR_dtheta * dtheta_dsigma
  df_dsigma = dR_dsigma * q + \
              R * dq_dsigma - dp_dsigma * sin_phi
  df_dsigma = df_dsigma / np.linalg.norm(df_dsigma)
  return df_dsigma

def dp_dsigma_fn():
  I = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).reshape(3,3)
  return I / 3.0

def dR_dtheta_fn(phi, theta):
  sqrt3 = np.sqrt(3)
  third = 1.0 / 3.0
  sin_phi = np.sin(phi * np.pi / 180)
  cos_phi = np.cos(phi * np.pi / 180)
  Theta = theta + third * np.pi
  sin_Theta = np.sin(Theta)
  cos_Theta = np.cos(Theta)
  R = 1 / sqrt3 * sin_Theta - third * cos_Theta * sin_phi
  dR_dtheta = 1 / sqrt3 * cos_Theta + third * sin_Theta * sin_phi
  return R, dR_dtheta

def dq_dsigma_fn(s, J2):
  sqrt3 = np.sqrt(3)
  dJ2_dsigma = dJ2_dsigma_fn(s)
  dq_dsigma = sqrt3 / (2.0 * np.sqrt(J2)) * dJ2_dsigma
  return dq_dsigma

def dJ2_dsigma_fn(s):
  return s

def dJ3_dsigma_fn(ss, J2):
  I = np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).reshape(3,3)
  dJ3_dsigma = ss - I * (2.0 / 3.0 * J2)
  return dJ3_dsigma
  
def dtheta_dsigma_fn(s, ss, J2, q, J3, theta):
  sin3theta = np.sin(3.0 * theta)
  dq_dsigma = dq_dsigma_fn(s, J2)
  dJ3_dsigma = dJ3_dsigma_fn(ss, J2)
  q_cub = q * q * q 
  dtheta_dsigma = -9.0 / (2.0 * sin3theta * q_cub) * \
                 (dJ3_dsigma - 3.0 * J3 * dq_dsigma / q)
  return dq_dsigma, dtheta_dsigma

def mohr_coulomb(c, phi, low, high, numcells):

    rad = np.linspace(0.01, high, numcells)
    grid = pv.CylinderStructured(radius = rad, height = (high - low),
                                 center = [0.5*(high+low), 0.5*(high+low), 0.5*(high+low)],
                                 direction = [1, 1, 1],
                                 theta_resolution = numcells,
                                 z_resolution = np.int(1.8*numcells))

    vtkGrid = pv.UnstructuredGrid(grid)
    stresses = vtkGrid.points

    sin_phi = np.sin(phi * np.pi / 180.0)
    cos_phi = np.cos(phi * np.pi / 180.0)
    sqrt3 = np.sqrt(3.0)

    yield_fn = []
    for stress in stresses:

      stress_mat = np.array([stress[0], 0.0, 0.0, 0.0, stress[1], 0.0, 0.0, 0.0, stress[2]]).reshape(3,3)

      # From my wikipedia version
      p, q, theta = pqtheta(stress_mat)
      Theta = theta + np.pi / 3.0
      Rmc = 1 / sqrt3 * np.sin(Theta) - 1 / 3 * np.cos(Theta) * sin_phi
      f = Rmc * q - p * sin_phi - c * cos_phi

      # normalize
      ind = np.argsort(stress)
      i1 = ind[2]
      i3 = ind[0]
      f = f / (abs(stress[i1]) + abs(stress[i3]) + 2.0 * c)

      yield_fn.append(f)

    vtkGrid["data"] = np.array(yield_fn)

    return vtkGrid

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
  
def plotData2():

  fac = 0.1
  lower_lim = -5.0e4
  upper_lim = 1.0e5 * fac
  ball_rad = 2000 * fac
  cyl_rad = 1000 * fac
  arrow_len = 1.0e5 * fac 

  c = 1.0e4
  phi = 30
  grid = mohr_coulomb(c, phi, lower_lim, upper_lim, 10)
  contours = grid.contour([0.0])

  snew   = np.array([-9596.73, -7.97093e-11, 3.28791e-10, -7.97093e-11, -14144.7, -9.62086e-10, 3.28791e-10, -9.62086e-10, -14745.9]).reshape(3,3)
  df_dsigma = df_dsigma_fn(snew, phi)
  snew_eig = np.diagonal(snew)
  snew_ball = pv.Sphere(radius=ball_rad, center = snew_eig)
  s_cyl, s_cone = drawArrow(snew_eig, np.diagonal(df_dsigma), arrow_len, cyl_rad)

  snew4 = np.array([-10798, -9.96366e-11, 4.10989e-10, -9.96366e-11, -15146.3, -1.20261e-09, 4.10989e-10, -1.20261e-09, -13677.5]).reshape(3,3)
  df_dsigma = df_dsigma_fn(snew4, phi)
  snew4_eig = np.diagonal(snew4)
  snew4_ball = pv.Sphere(radius=ball_rad, center = snew4_eig)
  s4_cyl, s4_cone = drawArrow(snew4_eig, np.diagonal(df_dsigma), arrow_len, cyl_rad)

  blue = np.array([12/256, 238/256, 246/256])
  grey = np.array([189/256, 189/256, 189/256])
  yellow = np.array([255/256, 247/256, 0/256])

  pv.set_plot_theme('document')
  p = pv.Plotter(notebook=False)
  #p.add_mesh(grid, opacity=0.4)
  p.add_mesh(contours, opacity=0.4, color = blue, show_scalar_bar=False)
  p.add_mesh(snew_ball, color="#ff7Eaa", show_edges=False)
  p.add_mesh(s_cyl, color="#ff7Eaa", show_edges=False)
  p.add_mesh(s_cone, color="#ff7Eaa", show_edges=False)
  p.add_mesh(snew4_ball, color="#ff7E00", show_edges=False)
  p.add_mesh(s4_cyl, color="#ff7E00", show_edges=False)
  p.add_mesh(s4_cone, color="#ff7E00", show_edges=False)
  p.add_axes()
  p.add_bounding_box()
  p.show_grid()
  p.show()

  return

plotData2()



