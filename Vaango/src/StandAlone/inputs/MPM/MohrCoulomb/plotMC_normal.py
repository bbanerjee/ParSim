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

  print("I1 = ", I1, "p = ", p, "J2 = ", J2, "q = ", q, "J3 = ", J3, "R^3 = ", r_cubed)

  R, dR_dtheta, Theta, sin_Theta, cos_Theta = dR_dtheta_fn(phi, theta)
  print("cos(3theta_c) = ", cos3theta, "theta_c = ", theta, "theta = ", Theta)
  print("cos(theta) = ", cos_Theta, "sin(theta) = ", sin_Theta, "R_phi = ", R)
  print("dR_phi_dtheta = ", dR_dtheta)

  dq_dsigma, dtheta_dsigma = dtheta_dsigma_fn(s, ss, J2, q, J3, theta)

  print("dq_dsigma = ", dq_dsigma)
  print("dtheta_dsigma = ", dtheta_dsigma)
  dp_dsigma = dp_dsigma_fn()
  dR_dsigma = dR_dtheta * dtheta_dsigma
  print("dR_phi_dsigma = ", dR_dsigma)
  print("dR_phi_dsigma * q = ", dR_dsigma * q)
  print("R_phi * dq_dsigma = ", R * dq_dsigma)
  print("dp_dsigma * sin_phi = ", dp_dsigma * sin_phi)
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
  return R, dR_dtheta, Theta, sin_Theta, cos_Theta

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
  print("sin(3theta_c) =", sin3theta, "q^3 = ", q_cub)
  print("dJ3_dsigma = ", dJ3_dsigma)
  fac_1 = -9.0 / (2.0 * sin3theta * q_cub) 
  print("fac_1 = ", fac_1)
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

def toMatrix(vec):
  return np.array([vec[0], vec[3], vec[4], vec[3], vec[1], vec[5], vec[4], vec[5], vec[2]]).reshape(3,3)
  
def plotData2():

  fac = 0.1
  lower_lim = -5.0e4
  upper_lim = 1.0e5 * fac
  ball_rad = 2000 * fac
  cyl_rad = 1000 * fac
  arrow_len = 1.0e5 * fac 

  c = 1.0e4
  phi = 30
  grid = mohr_coulomb(c, phi, lower_lim, upper_lim, 50)
  contours = grid.contour([0.0])

  #snew   = np.array([-9596.73, -7.97093e-11, 3.28791e-10, -7.97093e-11, -14144.7, -9.62086e-10, 3.28791e-10, -9.62086e-10, -14745.9]).reshape(3,3)
  #df_dsigma = df_dsigma_fn(snew, phi)
  #snew_eig = np.diagonal(snew)
  #snew_ball = pv.Sphere(radius=ball_rad, center = snew_eig)
  #s_cyl, s_cone = drawArrow(snew_eig, np.diagonal(df_dsigma), arrow_len, cyl_rad)

  #snew4 = np.array([-10798, -9.96366e-11, 4.10989e-10, -9.96366e-11, -15146.3, -1.20261e-09, 4.10989e-10, -1.20261e-09, -13677.5]).reshape(3,3)
  #df_dsigma = df_dsigma_fn(snew4, phi)
  #snew4_eig = np.diagonal(snew4)
  #snew4_ball = pv.Sphere(radius=ball_rad, center = snew4_eig)
  #s4_cyl, s4_cone = drawArrow(snew4_eig, np.diagonal(df_dsigma), arrow_len, cyl_rad)

  #df_dsigma_CPP = toMatrix(np.array([-0.633091, 0.412266, -0.65513, 2.59117e-15, 0.00604987, 1.53632e-15]))
  #snew5 = np.array([86603.9, 5.14732e-10, 3.77864e-12, 5.14732e-10, 294453, 3.04243e-10, 3.77864e-12, 3.04243e-10, 86603.9]).reshape(3,3)
  #df_dsigma = df_dsigma_fn(snew5, phi)
  #print("df_dsigma = ", df_dsigma)
  #print("df_dsigma_CPP = ", df_dsigma_CPP)
  #snew5_eig = np.diagonal(snew5)
  #snew5_ball = pv.Sphere(radius=ball_rad, center = snew5_eig)
  #s5_cyl, s5_cone = drawArrow(snew5_eig, np.diagonal(df_dsigma), arrow_len, cyl_rad)
  #s5_cyl_CPP, s5_cone_CPP = drawArrow(snew5_eig, np.diagonal(df_dsigma_CPP), arrow_len, cyl_rad)

  #df_dsigma_CPP = toMatrix(np.array([-0.0424937, -0.914076, 0.301887, 0.26744, 1.32117e-07, -9.68111e-08]))
  #snew6 = np.array([-9296.9, 1110.57, 0.00164243, 1110.57, -12916.2, -0.000385837, 0.00164243, -0.000385837, -5048.44]).reshape(3,3)
  #df_dsigma = df_dsigma_fn(snew6, phi)
  #print("df_dsigma = ", df_dsigma)
  #print("df_dsigma_CPP = ", df_dsigma_CPP)
  #snew6_eig = np.diagonal(snew6)
  #snew6_ball = pv.Sphere(radius=ball_rad, center = snew6_eig)
  #s6_cyl, s6_cone = drawArrow(snew6_eig, np.diagonal(df_dsigma), arrow_len, cyl_rad)
  #s6_cyl_CPP, s6_cone_CPP = drawArrow(snew6_eig, np.diagonal(df_dsigma_CPP), arrow_len, cyl_rad)
  
  #df_dsigma_CPP = toMatrix(np.array([-0.0723199058966196, -0.906982710164428, 0.326434205353643, 0.256111116996115, 1.53283933944362e-07, -9.23144656544059e-08]))
  #snew7 = np.array([-14433.5, 295.662, 0.00102008, 295.662, -15397.1, -9.40987e-05, 0.00102008, -9.40987e-05, -11800.7]).reshape(3,3)
  #df_dsigma = df_dsigma_fn(snew7, phi)
  #print("df_dsigma = ", df_dsigma)
  #print("df_dsigma_CPP = ", df_dsigma_CPP)
  #snew7_eig = np.diagonal(snew7)
  #snew7_ball = pv.Sphere(radius=ball_rad, center = snew7_eig)
  #s7_cyl, s7_cone = drawArrow(snew7_eig, np.diagonal(df_dsigma), arrow_len, cyl_rad)
  #s7_cyl_CPP, s7_cone_CPP = drawArrow(snew7_eig, np.diagonal(df_dsigma_CPP), arrow_len, cyl_rad)

  #df_dsigma_CPP = toMatrix(np.array([-0.0723199, -0.906983, 0.326434, 0.256111, 1.53284e-07, -9.23145e-08]))
  #snew8 = np.array([-14486.9, 48.3381, 0.00199985, 48.3381, -14644.4, 1.1731e-05, 0.00199985, 1.1731e-05, -9333.14]).reshape(3,3)
  #df_dsigma = df_dsigma_fn(snew8, phi)
  #print("df_dsigma = ", df_dsigma)
  #print("df_dsigma_CPP = ", df_dsigma_CPP)
  #snew8_eig = np.diagonal(snew8)
  #snew8_ball = pv.Sphere(radius=ball_rad, center = snew8_eig)
  #s8_cyl, s8_cone = drawArrow(snew8_eig, np.diagonal(df_dsigma), arrow_len, cyl_rad)
  #s8_cyl_CPP, s8_cone_CPP = drawArrow(snew8_eig, np.diagonal(df_dsigma_CPP), arrow_len, cyl_rad)

  #df_dsigma_CPP = toMatrix(np.array([-0.906982710165378, -0.0723199058949562, 0.326434205353318, -0.256111116993226, 4.801524242618e-07, 1.01684448941657e-07]))
  #stress = np.array([-15966.6, -319.478, 0.000947633, -319.478, -14925.5, 0.000132001, 0.000947633, 0.000132001, -13529.6]).reshape(3,3)
  #df_dsigma = df_dsigma_fn(stress, phi)
  #print("df_dsigma = ", df_dsigma)
  #print("df_dsigma_CPP = ", df_dsigma_CPP)
  #s_eig = np.diagonal(stress)
  #s_ball = pv.Sphere(radius=ball_rad, center = s_eig)
  #s_cyl, s_cone = drawArrow(s_eig, np.diagonal(df_dsigma), arrow_len, cyl_rad)
  #s_cyl_CPP, s_cone_CPP = drawArrow(s_eig, np.diagonal(df_dsigma_CPP), arrow_len, cyl_rad)

  df_dsigma_CPP = toMatrix(np.array([0.0953234866379436, -0.928974341477548,  0.170692592289366,  0.314299450697999, 2.74459914604108e-08, -1.15664865046287e-07]))
  stress = np.array([-15241.3, 311.429, 0.000328045, 311.429, -16256.2, -0.000110158, 0.000328045, -0.000110158, -14391.4]).reshape(3,3)
  df_dsigma = df_dsigma_fn(stress, phi)
  print("df_dsigma = ", df_dsigma)
  print("df_dsigma_CPP = ", df_dsigma_CPP)
  s_eig = np.diagonal(stress)
  s_ball = pv.Sphere(radius=ball_rad, center = s_eig)
  s_cyl, s_cone = drawArrow(s_eig, np.diagonal(df_dsigma), arrow_len, cyl_rad)
  s_cyl_CPP, s_cone_CPP = drawArrow(s_eig, np.diagonal(df_dsigma_CPP), arrow_len, cyl_rad)

  blue = np.array([12/256, 238/256, 246/256])
  grey = np.array([189/256, 189/256, 189/256])
  yellow = np.array([255/256, 247/256, 0/256])

  pv.set_plot_theme('document')
  p = pv.Plotter(notebook=False)
  #p.add_mesh(grid, opacity=0.4)
  p.add_mesh(contours, opacity=0.4, color = blue, show_scalar_bar=False)

  #p.add_mesh(snew_ball, color="#ff7Eaa", show_edges=False)
  #p.add_mesh(s_cyl, color="#ff7Eaa", show_edges=False)
  #p.add_mesh(s_cone, color="#ff7Eaa", show_edges=False)

  #p.add_mesh(snew4_ball, color="#ff7E00", show_edges=False)
  #p.add_mesh(s4_cyl, color="#ff7E00", show_edges=False)
  #p.add_mesh(s4_cone, color="#ff7E00", show_edges=False)

  #p.add_mesh(snew5_ball, color="#ff7E00", show_edges=False)
  #p.add_mesh(s5_cyl, color="#ff7E00", show_edges=False)
  #p.add_mesh(s5_cone, color="#ff7E00", show_edges=False)
  #p.add_mesh(s5_cyl_CPP, color="#dd7E00", show_edges=False)
  #p.add_mesh(s5_cone_CPP, color="#dd7E00", show_edges=False)

  #p.add_mesh(snew6_ball, color="#ff7E00", show_edges=False)
  #p.add_mesh(s6_cyl, color="#ff7E00", show_edges=False)
  #p.add_mesh(s6_cone, color="#ff7E00", show_edges=False)
  #p.add_mesh(s6_cyl_CPP, color="#dd7E00", show_edges=False)
  #p.add_mesh(s6_cone_CPP, color="#dd7E00", show_edges=False)

  #p.add_mesh(snew7_ball, color="#FF8243", show_edges=False)
  #p.add_mesh(s7_cyl, color="#FF8243", show_edges=False)
  #p.add_mesh(s7_cone, color="#FF8243", show_edges=False)
  #p.add_mesh(s7_cyl_CPP, color="#5E8C31	", show_edges=False)
  #p.add_mesh(s7_cone_CPP, color="#5E8C31", show_edges=False)

  #p.add_mesh(snew8_ball, color="#0A7E8C", show_edges=False)
  #p.add_mesh(s8_cyl, color="#0A7E8C", show_edges=False)
  #p.add_mesh(s8_cone, color="#0A7E8C", show_edges=False)
  #p.add_mesh(s8_cyl_CPP, color="#FFC40C", show_edges=False)
  #p.add_mesh(s8_cone_CPP, color="#FFC40C", show_edges=False)

  p.add_mesh(s_ball, color="#0A7E8C", show_edges=False)
  p.add_mesh(s_cyl, color="#0A7E8C", show_edges=False)
  p.add_mesh(s_cone, color="#0A7E8C", show_edges=False)
  p.add_mesh(s_cyl_CPP, color="#FFC40C", show_edges=False)
  p.add_mesh(s_cone_CPP, color="#FFC40C", show_edges=False)

  p.add_axes()
  p.add_bounding_box()
  p.show_grid()
  p.show()

  return

plotData2()



