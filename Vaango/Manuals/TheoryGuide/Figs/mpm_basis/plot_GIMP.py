import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def create_grid(low, high, numCells):
  xx = np.linspace(low, high, numCells+1)
  yy = np.linspace(low, high, numCells+1)
  x, y = np.meshgrid(xx, yy)
  z = x*0 + y*0
  node_coords = np.c_[x.reshape(-1), y.reshape(-1), z.reshape(-1)]
  return xx, yy, node_coords

def create_node_balls(points):
  balls = []
  for point in points:
    ball = pv.Sphere(radius=0.06, center = point)
    balls.append(ball)

  return balls

def create_particles(coords_g, h_g, l_p):
  x_g = coords_g[0]
  y_g = coords_g[1]
  x_p = np.array([x_g - h_g - l_p,
                  x_g - h_g,
                  x_g - h_g + l_p,
                  x_g - l_p,
                  x_g, 
                  x_g + l_p,
                  x_g + h_g - l_p,
                  x_g + h_g,
                  x_g + h_g + l_p])
  y_p = np.array([y_g - h_g - l_p,
                  y_g - h_g,
                  y_g - h_g + l_p,
                  y_g - l_p,
                  y_g, 
                  y_g + l_p,
                  y_g + h_g - l_p,
                  y_g + h_g,
                  y_g + h_g + l_p])
  balls = []
  pulses = []
  for x in x_p:
    for y in y_p:
      center = [x, y, 0]
      ball = pv.Sphere(radius=0.1, center = center)
      balls.append(ball)
      center = [x, y, 0.025]
      pulse = pv.Cube(center = center, x_length = 2*l_p,
                      y_length = 2*l_p, z_length = 0.05)
      pulses.append(pulse)
  
  return balls, pulses, x_p, y_p
  
def shape_function():
  R = np.array([-1, 1, 1, -1, -1, -1, 1, 1]).reshape(2,4)
  return R

def gimp_basis_single(nx, ny, node_coords, cur_node, h_g, l_p):
  corner_node_mm = [0]
  corner_node_mp = [nx - 1]
  corner_node_pp = [nx*ny - 1]
  corner_node_pm = [nx*ny - nx]
  edge_node_bot = list(range(1, nx - 1))
  edge_node_top = list(range(nx*ny - nx + 1, nx*ny - 1))
  edge_node_left = list(range(nx, nx*ny - 2*nx + 1, ny))
  edge_node_right = list(range(2*nx-1, nx*ny - nx, ny))
  boundary_nodes = []
  boundary_nodes.extend(corner_node_mm)
  boundary_nodes.extend(corner_node_mp)
  boundary_nodes.extend(corner_node_pp)
  boundary_nodes.extend(corner_node_pm)
  boundary_nodes.extend(edge_node_bot)
  boundary_nodes.extend(edge_node_top)
  boundary_nodes.extend(edge_node_left)
  boundary_nodes.extend(edge_node_right)

  surf = -1
  data = pv.PolyData()
  if cur_node not in boundary_nodes:
    x_g = node_coords[cur_node][0]
    y_g = node_coords[cur_node][1]
    x_support = np.linspace(x_g - h_g - l_p, x_g + h_g + l_p, 100)
    y_support = np.linspace(y_g - h_g - l_p, y_g + h_g + l_p, 100)
    vertices = compute_S_gimp(x_support, y_support, x_g, y_g, l_p, h_g)
    data.points = np.array(vertices) 
    surf = data.delaunay_2d();

  return surf

def compute_S_gimp(x_support, y_support, x_g, y_g, l_p, h_g):
  vertices = []
  S_gp = gimp_S_xy_gp(x_support, y_support, x_g, y_g, l_p, h_g)
  for ix, x in enumerate(x_support):
    for iy, y in enumerate(y_support):
       vertices.append([x, y, S_gp[ix, iy]])
  return np.array(vertices)

def gimp_S_xy_gp(x_p, y_p, x_g, y_g, l_p, h_g):
  S_x_gp = gimp_S_x_gp(x_p, x_g, l_p, h_g).reshape(x_p.shape[0],1)
  S_y_gp = gimp_S_x_gp(y_p, y_g, l_p, h_g).reshape(y_p.shape[0],1)
  S_xy_gp = np.matmul(S_x_gp, S_y_gp.transpose())
  return S_xy_gp

def gimp_S_x_gp(x_p, x_g, l_p, h_g):
  # Uintah implementation
  #
  #S_gp = [] 
  #for x_p_val in x_p:
  #  pg = (x_p_val - x_g) / h_g
  #  if (pg <= l_p):
  #    S_gp.append(1.0 - (pg * pg + l_p * l_p) / (2.0 * l_p))
  #  elif (pg > l_p and pg <= (1.0 - l_p)):
  #    S_gp.append(1.0 - pg)
  #  elif (pg > (1.0 - l_p)):
  #    S_gp.append((1.0 + l_p - pg) * (1.0 + l_p - pg) / (4.0 * l_p))
  #
  # Uintah manual
  #
  S_gp = [] 
  for x_p_val in x_p:
    pg = (x_p_val - x_g)
    if pg > -(h_g + l_p) and pg <= (-h_g + l_p):
      S_gp.append((h_g + l_p + pg) * (h_g + l_p + pg) / (4.0 * h_g * l_p))
    elif pg > (-h_g + l_p) and pg <= -l_p:
      S_gp.append(1 + pg / h_g)
    elif pg > -l_p and pg <= l_p:
      S_gp.append(1.0 - (pg * pg + l_p * l_p) / (2.0 * h_g * l_p))
    elif pg <= (h_g - l_p) and pg > l_p:
      S_gp.append(1 - pg / h_g)
    elif pg > (h_g - l_p) and pg <= (h_g + l_p):
      S_gp.append((h_g + l_p - pg) * (h_g + l_p - pg) / (4.0 * h_g * l_p))
    else:
      S_gp.append(0)

  #S_gp = list(map(lambda S : S * 2.0, S_gp))

  return np.array(S_gp) 

def linear_S_xy(x, y, x_g, y_g, l_p):
  S_x = linear_S_x(x, x_g, l_p).reshape(x.shape[0],1)
  S_y = linear_S_x(y, y_g, l_p).reshape(y.shape[0],1)
  S_xy = np.matmul(S_x, S_y.transpose())
  return S_xy

def linear_S_x(x, x_g, l_p):
  xi = map(lambda x_val : x_to_xi(x_val, x_g - l_p, x_g)  if (x_val < x_g) else \
                          x_to_xi(x_val, x_g, x_g+l_p), x)
  xi_g = map(lambda x_val : -1  if (x_val < x_g) else 1, x)
  xi = list(xi)
  xi_g = list(xi_g)
  S_x = linear_S_xi(xi, xi_g)
  return S_x

def x_to_xi(x, x_lo, x_hi):
  xi = 2.0 * (x - x_lo)/(x_hi - x_lo) - 1.0
  return xi

def linear_S_xi(xi, xi_g):
  S_xi = map(lambda xi_val, xi_g_val : \
               0.5 * (1 - xi_val * xi_g_val) \
               if (xi_val >= -1.0 or xi_val <= 1.0) else 0, xi, xi_g)
  return np.array(list(S_xi))

def compute_S_linear(x_support, y_support, x_g, y_g, l_p):
  vertices = []
  S_xy = linear_S_xy(x_support, y_support, x_g, y_g, l_p)
  for ix, x in enumerate(x_support):
    for iy, y in enumerate(y_support):
       vertices.append([x, y, S_xy[ix, iy]])
  return np.array(vertices)

def mpm_basis_single(nx, ny, node_coords, cur_node, h_g):
  corner_node_mm = [0]
  corner_node_mp = [nx - 1]
  corner_node_pp = [nx*ny - 1]
  corner_node_pm = [nx*ny - nx]
  edge_node_bot = list(range(1, nx - 1))
  edge_node_top = list(range(nx*ny - nx + 1, nx*ny - 1))
  edge_node_left = list(range(nx, nx*ny - 2*nx + 1, ny))
  edge_node_right = list(range(2*nx-1, nx*ny - nx, ny))
  boundary_nodes = []
  boundary_nodes.extend(corner_node_mm)
  boundary_nodes.extend(corner_node_mp)
  boundary_nodes.extend(corner_node_pp)
  boundary_nodes.extend(corner_node_pm)
  boundary_nodes.extend(edge_node_bot)
  boundary_nodes.extend(edge_node_top)
  boundary_nodes.extend(edge_node_left)
  boundary_nodes.extend(edge_node_right)

  surf = -1
  data = pv.PolyData()
  if cur_node not in boundary_nodes:
    x_g = node_coords[cur_node][0]
    y_g = node_coords[cur_node][1]
    x_support = np.linspace(x_g - h_g, x_g + h_g, 100)
    y_support = np.linspace(y_g - h_g, y_g + h_g, 100)
    vertices = compute_S_linear(x_support, y_support, x_g, y_g, h_g)
    data.points = np.array(vertices) 
    surf = data.delaunay_2d();

  return surf

if __name__ == "__main__":  

  low = 0.0
  high = 3.0
  nx = 3

  tab_cmap = plt.cm.get_cmap("autumn_r", nx*nx)

  pv.set_plot_theme('document')
  p = pv.Plotter(notebook=False)

  xx, yy, node_coords = create_grid(low, high, nx-1)
  balls = create_node_balls(node_coords)
  for ball in balls:
    p.add_mesh(ball, color="red", show_edges=False)

  h_g = (high - low)/(nx-1)
  l_p = 0.7 * h_g
  x_g = node_coords[4]
  particles, pulses, xpart, ypart = create_particles(x_g, h_g, l_p)

  for particle in particles:
    p.add_mesh(particle, color="green", show_edges=False)

  #for pulse in pulses:
  #  p.add_mesh(pulse, color=tab_cmap(0.4)[0:3], opacity = 0.2, show_edges=False)

  opacity = 0.8;
  color = 0;
  for node, coord in enumerate(node_coords):
    color = color + 1.0/(nx*nx)
    if (color >= 1.0):
      color = 0
    grid = mpm_basis_single(nx, nx, node_coords, node, h_g)
    p_grid = gimp_basis_single(nx, nx, node_coords, node, h_g, l_p)

    #if (grid != -1):
      #p.add_mesh(grid, color=tab_cmap(color)[0:3], opacity=opacity, 
    #  p.add_mesh(grid, color="blue", opacity=0.8, 
    #             specular=1.0, smooth_shading=False, show_edges=False)

    if (p_grid != -1):
      p.add_mesh(p_grid, color=tab_cmap(color)[0:3], opacity=1.0, 
                 specular=1.0, smooth_shading=False, show_edges=False)

  p.add_axes()
  p.add_bounding_box()
  #p.show_grid()
  p.show()

