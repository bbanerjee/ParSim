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
  sin3theta_s = r_cubed / (q * q * q)

  if (sin3theta_s < -1):
    sin3theta_s = -1
  if (sin3theta_s > 1):
    sin3theta_s = 1
  #print(sin3theta_s)
  theta_s = np.arcsin(sin3theta_s) / 3.0
  if (theta_s < -np.pi/6 or theta_s > np.pi/6):
    print(theta_s, stress)

  theta = np.pi/6 - theta_s
  return p, q, theta

def sigmbartheta(stress):

  I1, J2, J3 = I1J2J3(stress)

  sigma_m = I1 / 3.0
  sigma_bar = np.sqrt(J2)
  sin3theta = - 1.5 * np.sqrt(3) * J3 / (sigma_bar * sigma_bar * sigma_bar)
  if (sin3theta < -1):
    sin3theta = -1
  if (sin3theta > 1):
    sin3theta = 1
  #print(sin3theta)
  theta = np.arcsin(sin3theta) / 3.0

  return sigma_m, sigma_bar, theta

def mohr_coulomb(c, phi, low, high, numcells):

    rad = np.linspace(0.01, high, numcells)
    grid_abq = pv.CylinderStructured(radius = rad, height = (high - low),
                                 center = [0.5*(high+low), 0.5*(high+low), 0.5*(high+low)],
                                 direction = [1, 1, 1],
                                 theta_resolution = numcells,
                                 z_resolution = np.int(1.8*numcells))
    grid_mine = pv.CylinderStructured(radius = rad, height = (high - low),
                                 center = [0.5*(high+low), 0.5*(high+low), 0.5*(high+low)],
                                 direction = [1, 1, 1],
                                 theta_resolution = numcells,
                                 z_resolution = np.int(1.8*numcells))

    vtkGrid_abq = pv.UnstructuredGrid(grid_abq)
    vtkGrid_mine = pv.UnstructuredGrid(grid_mine)
    stresses = vtkGrid_abq.points

    sin_phi = np.sin(phi * np.pi / 180.0)
    cos_phi = np.cos(phi * np.pi / 180.0)
    tan_phi = np.tan(phi * np.pi / 180.0)
    sqrt3 = np.sqrt(3.0)

    yield_fn_abq = []
    yield_fn_mine = []
    for stress in stresses:

      stress_mat = np.array([stress[0], 0.0, 0.0, 0.0, stress[1], 0.0, 0.0, 0.0, stress[2]]).reshape(3,3)

      # From Abaqus manual
      p, q, theta = pqtheta(stress_mat)
      Theta = theta + np.pi / 3.0
      Rmc = 1 / (sqrt3 * cos_phi) * np.sin(Theta) + 1 / 3 * np.cos(Theta) * tan_phi
      f_abq = Rmc * q - p * tan_phi - c

      # From my wikipedia version
      #p, q, theta = pqtheta(stress_mat)
      #Theta = theta + np.pi / 3.0
      Rmc = 1 / sqrt3 * np.sin(Theta) - 1 / 3 * np.cos(Theta) * sin_phi
      f_mine = Rmc * q - p * sin_phi - c * cos_phi

      # From 1986 paper
      #sigm, sigbar, theta = sigmbartheta(stress_mat)
      #Theta = theta
      #Rmc = np.cos(Theta) - 1/sqrt3 * np.cos(Theta) * sin_phi
      #f = Rmc * sigbar - sigm * sin_phi - c * cos_phi

      # normalize
      ind = np.argsort(stress)
      i1 = ind[2]
      i3 = ind[0]
      f_abq = f_abq / (abs(stress[i1]) + abs(stress[i3]) + 2.0 * c)
      f_mine = f_mine / (abs(stress[i1]) + abs(stress[i3]) + 2.0 * c)

      yield_fn_abq.append(f_abq)
      yield_fn_mine.append(f_mine)

    vtkGrid_abq["data"] = np.array(yield_fn_abq)
    vtkGrid_mine["data"] = np.array(yield_fn_mine)

    return vtkGrid_abq, vtkGrid_mine

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

  fac = 1
  lower_lim = -5.0e4 * fac
  upper_lim = 1.0e5 * fac
  ball_rad = 2000 * fac
  cyl_rad = 1000 * fac
  arrow_len = 1.0e5 * fac 

  c = 1.0e4
  phi = 30
  grid_abq, grid_mine = mohr_coulomb(c, phi, lower_lim, upper_lim, 100)
  contours_abq = grid_abq.contour([0.0])
  contours_mine = grid_mine.contour([0.0])

  blue = np.array([12/256, 238/256, 246/256])
  grey = np.array([189/256, 189/256, 189/256])
  yellow = np.array([255/256, 247/256, 0/256])

  pv.set_plot_theme('document')
  p = pv.Plotter(notebook=False)
  #p.add_mesh(grid, opacity=0.4)
  p.add_mesh(contours_abq, opacity=0.7, color = "red", show_scalar_bar=False)
  p.add_mesh(contours_mine, opacity=0.4, color = blue, show_scalar_bar=False)
  p.add_axes()
  p.add_bounding_box()
  p.show_grid()
  p.show()

  return

plotData2()



