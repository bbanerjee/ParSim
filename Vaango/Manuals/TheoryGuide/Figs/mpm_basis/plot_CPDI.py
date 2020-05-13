import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def create_grid(low, high, numCells):
  xx = np.linspace(low, high, numCells+1)
  yy = np.linspace(low, high, numCells+1)
  x, y = np.meshgrid(xx, yy)
  z = np.ones(x.shape)
  node_coords = np.c_[x.reshape(-1), y.reshape(-1), z.reshape(-1)]
  return xx, yy, node_coords

def create_grid_lines(xx, yy):
  nx = xx.shape[0]
  ny = yy.shape[0]
  dx = xx[1] - xx[0]
  dy = yy[1] - yy[0]
  grid_lines = []
  for x in xx:
    cyl_bot = np.array([x, yy[0] - dy, 1.0])
    cyl_top = np.array([x, yy[ny-1] + dy, 1.0])
    cyl_cen = 0.5 * (cyl_bot + cyl_top)
    direction = cyl_top - cyl_bot
    height = np.linalg.norm(direction)
    direction = direction / height
    radius = 0.01
    grid_cyl = pv.Cylinder(cyl_cen, direction = direction, radius = radius, 
                          height = height)
    grid_lines.append(grid_cyl)

  for y in yy:
    cyl_bot = np.array([xx[0] - dx, y, 1.0])
    cyl_top = np.array([xx[nx-1] + dx, y, 1.0])
    cyl_cen = 0.5 * (cyl_bot + cyl_top)
    direction = cyl_top - cyl_bot
    height = np.linalg.norm(direction)
    direction = direction / height
    radius = 0.01
    grid_cyl = pv.Cylinder(cyl_cen, direction = direction, radius = radius, 
                          height = height)
    grid_lines.append(grid_cyl)

  return grid_lines
  
def create_node_balls(points):
  balls = []
  for point in points:
    ball = pv.Sphere(radius=0.06, center = point)
    balls.append(ball)

  return balls

def create_particles(low, high, coords_g, h_g, l_p):
  parts = np.array([9.43699879e-04, 1, 1, 3940752754737155, 1.17336686e+00, -3.05365323e-02, -5.58774293e-02, -1.41675709e-02, 1.12901402e+00, -5.78028055e-02, -4.94209329e-01, -4.86743641e-01, 7.06774612e-01, 4.94356845e-01, 4.94282113e-01, -4.38364657e-02,
                    9.43699879e-04, 2, 1, 3940757049704450, 1.17232565e+00, 3.42853785e-02, -5.57971092e-02, 1.86067603e-02, 1.13204601e+00, 5.80542914e-02, -4.92291200e-01, 4.79577565e-01, 7.06523735e-01, 4.94375296e-01, 5.05571443e-01, -4.38759311e-02,
                    9.43699879e-04, 2, 1, 3940757049704451, 1.15132661e+00, 3.48080696e-02, -4.29509574e-02, 5.19450276e-02, 2.36043360e+00, 2.75139040e-01, -4.87537509e-01, 1.19416547e+00, 6.91743232e-01, 4.94776970e-01, 5.23026948e-01, -3.55605706e-02,
                    9.43699879e-04, 3, 1, 3940752754802689, 1.15872247e+00, -1.52954316e-02, 6.42199131e-02, -5.34594779e-03, 1.15273671e+00, -6.01853570e-02, 5.23718188e-01, -5.06278797e-01, 7.12086924e-01, 5.05991277e-01, 4.94182737e-01, -4.36540845e-02,
                    9.43699879e-04, 3, 1, 3940752754802691, 2.38035320e+00, 3.59929705e-03, 2.73277007e-01, -2.41068750e-03, 1.15380733e+00, -5.06575537e-02, 1.23549735e+00, -5.37361644e-01, 6.92272195e-01, 5.23683767e-01, 4.94186033e-01, -3.49000328e-02,
                    9.43699879e-04, 4, 1, 3940757049769984, 1.15731141e+00, 1.66114819e-02, 6.42297460e-02, 5.99177653e-03, 1.15542969e+00, 6.06312366e-02, 5.23514190e-01, 5.00859949e-01, 7.11833044e-01, 5.05998020e-01, 5.05696379e-01, -4.36842927e-02,
                    9.43699879e-04, 4, 1, 3940757049769985, 1.13084485e+00, 2.05308081e-02, 5.48862791e-02, 1.25800255e-02, 2.37751841e+00, 2.74840084e-01, 5.60758344e-01, 1.20832842e+00, 6.94343862e-01, 5.06152320e-01, 5.23368253e-01, -3.51515053e-02,
                    9.43699879e-04, 4, 1, 3940757049769986, 2.37983813e+00, -4.17582891e-03, 2.73601316e-01, 3.12704960e-03, 1.15751862e+00, 5.11926547e-02, 1.23359860e+00, 5.30612240e-01, 6.90917997e-01, 5.23679682e-01, 5.05699093e-01, -3.49366480e-02,
                    9.43699879e-04, 4, 1, 3940757049769987, 2.23115603e+00, -6.72744852e-02, 2.19768609e-01, -6.12256826e-02, 2.25162131e+00, 2.22127726e-01, 1.22910403e+00, 1.20384610e+00, 6.45629079e-01, 5.22808849e-01, 5.22589676e-01, -2.77843328e-02])
  parts = parts.reshape(9, 16)
  particles = []
  shapes = []
  xps = []
  xp1s = []
  xp2s = []
  xp3s = []
  xp4s = []
  for part in parts:
    t = (0.515 - part[13])/(0.515 - 0.485) + 0.2
    x = (1 - t) * low + t * high
    t = (0.515 - part[14])/(0.515 - 0.485) + 0.2
    y = (1 - t) * low + t * high
    center = [x, y, -0.1]
    xps.append(np.array(center))
    F = part[4:13].reshape(3,3)
    xp1, xp2, xp3, xp4, surf = deform_shape(x, y, F, 0.5*h_g)
    xp1s.append(xp1)
    xp2s.append(xp2)
    xp3s.append(xp3)
    xp4s.append(xp4)
    shapes.append(surf)
    particle = pv.Sphere(radius=0.06, center = center)
    particles.append(particle)
  
  return particles, shapes, xps, xp1s, xp2s, xp3s, xp4s

def deform_shape(x, y, F, l):
  data = pv.PolyData()
  r1_init = np.array([0.65*l, 0, 0])
  r2_init = np.array([0, 0.65*l, 0])
  F[0,2] = 0
  F[1,2] = 0
  F[2,0] = 0
  F[2,1] = 0
  F[2,2] = 1
  R = np.column_stack((r1_init, r2_init))
  R_deformed = np.matmul(F, R)
  r1 = R_deformed[:,0]
  r2 = R_deformed[:,1]
  xp = np.array([x, y, -0.1]).transpose()
  xp_1 = xp - (r1+r2)
  xp_2 = xp + (r1-r2)
  xp_3 = xp + (r1+r2)
  xp_4 = xp - (r1-r2)
  vertices = [xp_1, xp_2, xp_3, xp_4]
  data.points = np.array(vertices) 
  surf = data.delaunay_2d();
  surf.extrude([0, 0, 1])
  return xp_1, xp_2, xp_3, xp_4, surf

def cpdi_S_xy(x, xp1, xp2, xp3, xp4):
  return 0 
  
def linear_xieta_to_xy(xi, eta, xp1, xp2, xp3, xp4):
  N1 = 1.0 / 4.0 * (1 - xi)*(1 - eta)
  N2 = 1.0 / 4.0 * (1 + xi)*(1 - eta)
  N3 = 1.0 / 4.0 * (1 + xi)*(1 + eta)
  N4 = 1.0 / 4.0 * (1 - xi)*(1 + eta)
  x = N1 * xp1 + N2 * xp2 + N3 * xp3 + N4 * xp4
  return x

def linear_xy_to_xieta(x, xp1, xp2, xp3, xp4):
  a0 = xp1[0] + xp2[0] + xp3[0] + xp4[0]
  a1 = xp1[1] + xp2[1] + xp3[1] + xp4[1]
  b0 = -xp1[0] + xp2[0] + xp3[0] - xp4[0]
  b1 = -xp1[1] + xp2[1] + xp3[1] - xp4[1]
  c0 = -xp1[0] - xp2[0] + xp3[0] + xp4[0]
  c1 = -xp1[1] - xp2[1] + xp3[1] + xp4[1]
  d0 = xp1[0] - xp2[0] + xp3[0] - xp4[0]
  d1 = xp1[1] - xp2[1] + xp3[1] - xp4[1]

  xieta = np.array([0, 0]).reshape(2,1)
  while True:
    f = np.array([a0 + b0*xieta[0] + c0*xieta[1] + d0*xieta[0]*xieta[1] - 4.0 * x[0],
                  a1 + b1*xieta[0] + c1*xieta[1] + d1*xieta[0]*xieta[1] - 4.0 * x[1]]).reshape(2,1)
    df = np.array([b0 + d0*xieta[1], c0 + d0*xieta[0],
                   b1 + d1*xieta[1], c1 + d1*xieta[0]]).reshape(2,2)
    dxieta = np.linalg.solve(df, f)
    xieta = xieta - dxieta
    if (np.linalg.norm(dxieta) < 1.0e-6):
      break;

  return xieta.flatten()

def linear_S_xieta(xi, eta, xieta_g):
  S_xieta = 1.0 / 4.0 * (1 + xi * xieta_g[0]) * (1 + eta * xieta_g[1])
  if (xi < -1.0 or xi > 1.0 or eta < -1.0 or eta > 1.0):
    S_xieta = 0.0
  return S_xieta

def compute_S_linear_xieta(xieta_g, xp1, xp2, xp3, xp4):
  xi_support = np.linspace(-1, 1, 50)
  eta_support = np.linspace(-1, 1, 50)
  vertices = []
  for xi in xi_support:
    for eta in eta_support:
      S_xieta = linear_S_xieta(xi, eta, xieta_g)
      xp = linear_xieta_to_xy(xi, eta, xp1, xp2, xp3, xp4) 
      vertices.append([xp[0], xp[1], S_xieta])
  return np.array(vertices)

def cpdi_Stilde_g(xieta, xieta_g, xg1, xg2, xg3, xg4, xp1, xp2, xp3, xp4):
  xieta_p_p1 = np.array([-1,-1])
  xieta_p_p2 = np.array([1,-1])
  xieta_p_p3 = np.array([1,1])
  xieta_p_p4 = np.array([-1,1])
  N_xieta_p_p1 = linear_S_xieta(xieta[0], xieta[1], xieta_p_p1)
  N_xieta_p_p2 = linear_S_xieta(xieta[0], xieta[1], xieta_p_p2)
  N_xieta_p_p3 = linear_S_xieta(xieta[0], xieta[1], xieta_p_p3)
  N_xieta_p_p4 = linear_S_xieta(xieta[0], xieta[1], xieta_p_p4)
  xieta_g_p1 = linear_xy_to_xieta(xp1, xg1, xg2, xg3, xg4)
  xieta_g_p2 = linear_xy_to_xieta(xp2, xg1, xg2, xg3, xg4)
  xieta_g_p3 = linear_xy_to_xieta(xp3, xg1, xg2, xg3, xg4)
  xieta_g_p4 = linear_xy_to_xieta(xp4, xg1, xg2, xg3, xg4)
  S_xieta_g_p1 = linear_S_xieta(xieta_g_p1[0], xieta_g_p1[1], xieta_g)
  S_xieta_g_p2 = linear_S_xieta(xieta_g_p2[0], xieta_g_p2[1], xieta_g)
  S_xieta_g_p3 = linear_S_xieta(xieta_g_p3[0], xieta_g_p3[1], xieta_g)
  S_xieta_g_p4 = linear_S_xieta(xieta_g_p4[0], xieta_g_p4[1], xieta_g)
  Stilde_g = (N_xieta_p_p1 * S_xieta_g_p1 + N_xieta_p_p2 * S_xieta_g_p2 + 
              N_xieta_p_p3 * S_xieta_g_p3 + N_xieta_p_p4 * S_xieta_g_p4)
  #print("x_g = ", xg1, "xp1 = ", xp1, "xp2 = ", xp2, "xp3 = ", xp3, "xp4 = ", xp4)
  #print("N_p1 = ", N_xieta_p_p1, "N_p2 = ", N_xieta_p_p2, "N_p3 = ", N_xieta_p_p3, "N_p4 = ", N_xieta_p_p4)
  #print("xi_g1 = ", xieta_g_p1, "xi_g2 = ", xieta_g_p2, "xi_g3 = ", xieta_g_p3, "xi_g4 = ", xieta_g_p4)
  #print("S_p1 = ", S_xieta_g_p1, "S_p2 = ", S_xieta_g_p2, "S_p3 = ", S_xieta_g_p3, "S_p4 = ", S_xieta_g_p4)
  #print("Stilde_g = ", Stilde_g)
  return Stilde_g

def compute_Stilde_cpdi_xieta(xieta_g, xg1, xg2, xg3, xg4, xp1, xp2, xp3, xp4):
  xi_support = np.linspace(-1, 1, 5)
  eta_support = np.linspace(-1, 1, 5)
  vertices = []
  for xi in xi_support:
    for eta in eta_support:
      xieta = np.array([xi, eta])
      Stilde_g = cpdi_Stilde_g(xieta, xieta_g, xg1, xg2, xg3, xg4, xp1, xp2, xp3, xp4)
      xg = linear_xieta_to_xy(xi, eta, xg1, xg2, xg3, xg4) 
      coord = [xg[0], xg[1], Stilde_g]
      #print(coord)
      vertices.append(coord)
  vertices = np.array(vertices)
  return vertices

def compute_Sg_xp(x_g, h_g, x_p): 
  Sg_p = 0
  if (x_p[0] - x_g[0] >= 0.0 and x_p[1] - x_g[1] >= 0.0): 
    xg1 = x_g
    xg2 = x_g + np.array([h_g, 0])
    xg3 = x_g + np.array([h_g, h_g])
    xg4 = x_g + np.array([0, h_g])
    xieta_g = np.array([-1,-1])
    xieta_g_p = linear_xy_to_xieta(x_p, xg1, xg2, xg3, xg4)
    Sg_p = linear_S_xieta(xieta_g_p[0], xieta_g_p[1], xieta_g)
  elif (x_p[0] - x_g[0] < 0.0 and  x_p[1] - x_g[1] >= 0.0): 
    xg2 = x_g
    xg3 = x_g + np.array([0, h_g])
    xg4 = x_g + np.array([-h_g, h_g])
    xg1 = x_g + np.array([-h_g, 0])
    xieta_g = np.array([1,-1])
    xieta_g_p = linear_xy_to_xieta(x_p, xg1, xg2, xg3, xg4)
    Sg_p = linear_S_xieta(xieta_g_p[0], xieta_g_p[1], xieta_g)
  elif (x_p[0] - x_g[0] < 0.0 and x_p[1] - x_g[1] < 0.0): 
    xg3 = x_g
    xg4 = x_g + np.array([-h_g, 0])
    xg1 = x_g + np.array([-h_g, -h_g])
    xg2 = x_g + np.array([0, -h_g])
    xieta_g = np.array([1,1])
    xieta_g_p = linear_xy_to_xieta(x_p, xg1, xg2, xg3, xg4)
    Sg_p = linear_S_xieta(xieta_g_p[0], xieta_g_p[1], xieta_g)
  elif (x_p[0] - x_g[0] >= 0.0 and x_p[1] - x_g[1] < 0.0): 
    xg4 = x_g
    xg1 = x_g + np.array([0, -h_g])
    xg2 = x_g + np.array([h_g, -h_g])
    xg3 = x_g + np.array([h_g, 0])
    xieta_g = np.array([-1,1])
    xieta_g_p = linear_xy_to_xieta(x_p, xg1, xg2, xg3, xg4)
    Sg_p = linear_S_xieta(xieta_g_p[0], xieta_g_p[1], xieta_g)

  return Sg_p

def compute_Np_x(xieta, xp1, xp2, xp3, xp4): 
  xieta_p1 = np.array([-1,-1])
  xieta_p2 = np.array([1,-1])
  xieta_p3 = np.array([1,1])
  xieta_p4 = np.array([-1,1])
  N1 = linear_S_xieta(xieta[0], xieta[1], xieta_p1)
  N2 = linear_S_xieta(xieta[0], xieta[1], xieta_p2)
  N3 = linear_S_xieta(xieta[0], xieta[1], xieta_p3)
  N4 = linear_S_xieta(xieta[0], xieta[1], xieta_p4)
  return N1, N2, N3, N4

def cpdi_basis_grid(x_g, h_g, xp1, xp2, xp3, xp4):
  Sg_p1 = compute_Sg_xp(x_g, h_g, xp1) 
  Sg_p2 = compute_Sg_xp(x_g, h_g, xp2) 
  Sg_p3 = compute_Sg_xp(x_g, h_g, xp3) 
  Sg_p4 = compute_Sg_xp(x_g, h_g, xp4) 
  xis = np.linspace(-1, 1, 40)
  etas = np.linspace(-1, 1, 40)
  vertices = []
  for xi in xis:
    for eta in etas:
      xieta = np.array([xi, eta])
      N1, N2, N3, N4 = compute_Np_x(xieta, xp1, xp2, xp3, xp4) 
      Stilde_g = (N1 * Sg_p1 + N2 * Sg_p2 + N3 * Sg_p3 + N4 * Sg_p4)
      xp = linear_xieta_to_xy(xi, eta, xp1, xp2, xp3, xp4) 
      coord = [xp[0], xp[1], Stilde_g]
      #print(coord)
      vertices.append(coord)
  data = pv.PolyData(np.array(vertices))
  surf = data.delaunay_2d()
  return surf

def mpm_basis_particle(xp1, xp2, xp3, xp4):
  xieta_g1 = np.array([-1,-1])
  xieta_g2 = np.array([1,-1])
  xieta_g3 = np.array([1,1])
  xieta_g4 = np.array([-1,1])
  vertices_g1 = compute_S_linear_xieta(xieta_g1, xp1, xp2, xp3, xp4)
  vertices_g2 = compute_S_linear_xieta(xieta_g2, xp1, xp2, xp3, xp4)
  vertices_g3 = compute_S_linear_xieta(xieta_g3, xp1, xp2, xp3, xp4)
  vertices_g4 = compute_S_linear_xieta(xieta_g4, xp1, xp2, xp3, xp4)
  data_g1 = pv.PolyData(vertices_g1)
  data_g2 = pv.PolyData(vertices_g2)
  data_g3 = pv.PolyData(vertices_g3)
  data_g4 = pv.PolyData(vertices_g4)
  surf_g1 = data_g1.delaunay_2d()
  surf_g2 = data_g2.delaunay_2d()
  surf_g3 = data_g3.delaunay_2d()
  surf_g4 = data_g4.delaunay_2d()
  return surf_g1, surf_g2, surf_g3, surf_g4

def mpm_basis_grid(x_g, h_g):
  xg1 = x_g
  xg2 = x_g + np.array([h_g, 0])
  xg3 = x_g + np.array([h_g, h_g])
  xg4 = x_g + np.array([0, h_g])
  xieta_g = np.array([-1,-1])
  vertices_1 =  compute_S_linear_xieta(xieta_g, xg1, xg2, xg3, xg4)
  xg2 = x_g
  xg3 = x_g + np.array([0, h_g])
  xg4 = x_g + np.array([-h_g, h_g])
  xg1 = x_g + np.array([-h_g, 0])
  xieta_g = np.array([1,-1])
  vertices_2 =  compute_S_linear_xieta(xieta_g, xg1, xg2, xg3, xg4)
  xg3 = x_g
  xg4 = x_g + np.array([-h_g, 0])
  xg1 = x_g + np.array([-h_g, -h_g])
  xg2 = x_g + np.array([0, -h_g])
  xieta_g = np.array([1,1])
  vertices_3 =  compute_S_linear_xieta(xieta_g, xg1, xg2, xg3, xg4)
  xg4 = x_g
  xg1 = x_g + np.array([0, -h_g])
  xg2 = x_g + np.array([h_g, -h_g])
  xg3 = x_g + np.array([h_g, 0])
  xieta_g = np.array([-1,1])
  vertices_4 =  compute_S_linear_xieta(xieta_g, xg1, xg2, xg3, xg4)
  data_g1 = pv.PolyData(vertices_1)
  data_g2 = pv.PolyData(vertices_2)
  data_g3 = pv.PolyData(vertices_3)
  data_g4 = pv.PolyData(vertices_4)
  surf_g1 = data_g1.delaunay_2d()
  surf_g2 = data_g2.delaunay_2d()
  surf_g3 = data_g3.delaunay_2d()
  surf_g4 = data_g4.delaunay_2d()
  return surf_g1, surf_g2, surf_g3, surf_g4


if __name__ == "__main__":  

  low = 0.0
  high = 3.0
  nx = 3

  pv.set_plot_theme('document')
  p = pv.Plotter(notebook=False)

  xx, yy, node_coords = create_grid(low, high, nx-1)
  balls = create_node_balls(node_coords)
  grid_lines = create_grid_lines(xx, yy)
  #for ball in balls:
  #  p.add_mesh(ball, color="red", show_edges=False)
  for line in grid_lines:
    p.add_mesh(line, color="black", show_edges=False)

  tab_cmap = plt.cm.get_cmap("tab20", nx*nx)
  colors =  plt.cm.tab20(np.arange(10).astype(int))
  #print(colors)
  #plt.scatter(np.arange(10),np.ones(10), c=colors, s=180)
  #plt.show()

  h_g = (high - low)/(nx-1)
  l_p = 0.7 * h_g
  x_g = node_coords[4]

  particles, pulses, xps, xp1s, xp2s, xp3s, xp4s = create_particles(low, high, x_g, h_g, l_p)

  for particle in particles:
    p.add_mesh(particle, color="green", show_edges=False)

  color = 0
  for pulse in pulses:
    color = color + 1
    if (color >= 10):
      color = 0
    p.add_mesh(pulse, color=colors[color][0:3], opacity = 1.0, show_edges=False)

  opacity = 1.0
  color = 0
  for coord in node_coords:
    #coord = node_coords[4]
    color = color + 1
    if (color >= 10):
      color = 0
    x_g = np.array([coord[0], coord[1]])
    node = pv.Sphere(radius=0.06, center = [x_g[0], x_g[1], 1.0])
    p.add_mesh(node, color=colors[color][0:3], show_edges=False)
    for xp1, xp2, xp3, xp4 in zip(xp1s, xp2s, xp3s, xp4s):
      #xp1 = xp1s[8]
      #p2 = xp2s[8]
      #p3 = xp3s[8]
      #p4 = xp4s[8]
      d_p1_g = np.abs(x_g - xp1[0:2])
      d_p2_g = np.abs(x_g - xp2[0:2])
      d_p3_g = np.abs(x_g - xp3[0:2])
      d_p4_g = np.abs(x_g - xp4[0:2])
      if (d_p1_g[0] <= h_g and d_p1_g[1] <= h_g):
        particle = pv.Sphere(radius=0.06, center = [xp1[0], xp1[1], 0.0])
        p.add_mesh(particle, color="yellow", show_edges=False)
      if (d_p2_g[0] <= h_g and d_p2_g[1] <= h_g):
        particle = pv.Sphere(radius=0.06, center = [xp2[0], xp2[1], 0.0])
        p.add_mesh(particle, color="yellow", show_edges=False)
      if (d_p3_g[0] <= h_g and d_p3_g[1] <= h_g):
        particle = pv.Sphere(radius=0.06, center = [xp3[0], xp3[1], 0.0])
        p.add_mesh(particle, color="yellow", show_edges=False)
      if (d_p4_g[0] <= h_g and d_p4_g[1] <= h_g):
        particle = pv.Sphere(radius=0.06, center = [xp4[0], xp4[1], 0.0])
        p.add_mesh(particle, color="yellow", show_edges=False)

      surf_g = cpdi_basis_grid(x_g, h_g, xp1, xp2, xp3, xp4)
      p.add_mesh(surf_g, color=colors[color][0:3], opacity=opacity, 
                 specular=1.0, smooth_shading=False, show_edges=False)

  #opacity = 1.0
  #for pulse, xp1, xp2, xp3, xp4 in zip(pulses, xp1s, xp2s, xp3s, xp4s):
  #  color = color + 0.1
  #  if (color >= 1.0):
  #    color = 0
  #  p.add_mesh(pulse, color=tab_cmap(color)[0:3], opacity = 0.3, show_edges=False)
  #  surf_g1, surf_g2, surf_g3, surf_g4 = mpm_basis_particle(xp1, xp2, xp3, xp4)
  #  p.add_mesh(surf_g1, color=tab_cmap(color)[0:3], opacity=opacity, 
  #             specular=1.0, smooth_shading=False, show_edges=False)
  #  p.add_mesh(surf_g2, color=tab_cmap(color)[0:3], opacity=opacity, 
  #             specular=1.0, smooth_shading=False, show_edges=False)
  #  p.add_mesh(surf_g3, color=tab_cmap(color)[0:3], opacity=opacity, 
  #             specular=1.0, smooth_shading=False, show_edges=False)
  #  p.add_mesh(surf_g4, color=tab_cmap(color)[0:3], opacity=opacity, 
  #             specular=1.0, smooth_shading=False, show_edges=False)

  #opacity = 1.0
  #for coord in node_coords:
  #  color = color + 0.1
  #  if (color >= 1.0):
  #    color = 0
  #  x_g = np.array([coord[0], coord[1]])
  #  surf_g1, surf_g2, surf_g3, surf_g4 = mpm_basis_grid(x_g, h_g)
  #  p.add_mesh(surf_g1, color=tab_cmap(color)[0:3], opacity=opacity, 
  #             specular=1.0, smooth_shading=False, show_edges=False)
  #  p.add_mesh(surf_g2, color=tab_cmap(color)[0:3], opacity=opacity, 
  #             specular=1.0, smooth_shading=False, show_edges=False)
  #  p.add_mesh(surf_g3, color=tab_cmap(color)[0:3], opacity=opacity, 
  #             specular=1.0, smooth_shading=False, show_edges=False)
  #  p.add_mesh(surf_g4, color=tab_cmap(color)[0:3], opacity=opacity, 
  #             specular=1.0, smooth_shading=False, show_edges=False)

  p.add_axes()
  #p.add_bounding_box()
  #p.show_grid()
  p.show()

