import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def create_grid(low, high, numCells):
  xx = np.linspace(low, high, numCells+1)
  yy = np.linspace(low, high, numCells+1)
  x, y = np.meshgrid(xx, yy)
  z = x*0 + y*0
  points = np.c_[x.reshape(-1), y.reshape(-1), z.reshape(-1)]
  return x, y, points

def create_node_balls(points):
  balls = []
  for point in points:
    ball = pv.Sphere(radius=0.03, center = point)
    balls.append(ball)

  return balls

def create_particle_balls(low, high, numParticles):
  np.random.seed(12345)
  numCells = np.int(np.sqrt(numParticles))
  xx = np.linspace(low, high, numCells+1)
  yy = np.linspace(low, high, numCells+1)
  x, y = np.meshgrid(xx, yy)
  z = x*0 + y*0
  points = np.c_[x.reshape(-1), y.reshape(-1), z.reshape(-1)]
  particles = []
  for point in points:
    x_move = np.random.uniform(0, 0.1*high, size = 1)
    y_move = np.random.uniform(0, 0.1*high, size = 1)
    center = [point[0] + x_move[0], point[1] + y_move[0], 0]
    particle = pv.Sphere(radius=0.03, center = center)
    particles.append(particle)
  
  return particles, xx, yy
  
def shape_function():
  R = np.array([-1, 1, 1, -1, -1, -1, 1, 1]).reshape(2,4)
  return R

def linear_S_xy(x, y, x_g, y_g, l_g):
  S_x = linear_S_x(x, x_g, l_g).reshape(x.shape[0],1)
  S_y = linear_S_x(y, y_g, l_g).reshape(y.shape[0],1)
  S_xy = np.matmul(S_x, S_y.transpose())
  return S_xy

def linear_S_x(x, x_g, l_g):
  xi = map(lambda x_val : x_to_xi(x_val, x_g - l_g, x_g)  if (x_val < x_g) else \
                          x_to_xi(x_val, x_g, x_g+l_g), x)
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

def compute_S(x_support, y_support, x_g, y_g, l_g):
  vertices = []
  S_xy = linear_S_xy(x_support, y_support, x_g, y_g, l_g)
  for ix, x in enumerate(x_support):
    for iy, y in enumerate(y_support):
       vertices.append([x, y, S_xy[ix, iy]])
  return np.array(vertices)

def mpm_basis_actual(nx, ny, node_coords, cur_node):
  l_g = (node_coords[1][0] - node_coords[0][0])
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
  #print(boundary_nodes)

  vertices = []
  for node, coord in enumerate(node_coords):
    #print("cur_node = ", cur_node, "node = ", node)
    if cur_node == node:
      x_g = coord[0]
      y_g = coord[1]
      if node in corner_node_mm:
        x_support = np.linspace(x_g, x_g + l_g, 50)
        y_support = np.linspace(y_g, y_g + l_g, 50)
        vertices = compute_S(x_support, y_support, x_g, y_g, l_g)
      if cur_node in corner_node_mp:
        x_support = np.linspace(x_g - l_g, x_g , 50)
        y_support = np.linspace(y_g, y_g + l_g, 50)
        vertices = compute_S(x_support, y_support, x_g, y_g, l_g)
      if cur_node in corner_node_pm:
        x_support = np.linspace(x_g, x_g + l_g, 50)
        y_support = np.linspace(y_g - l_g, y_g, 50)
        vertices = compute_S(x_support, y_support, x_g, y_g, l_g)
      if cur_node in corner_node_pp:
        x_support = np.linspace(x_g - l_g, x_g, 50)
        y_support = np.linspace(y_g - l_g, y_g, 50)
        vertices = compute_S(x_support, y_support, x_g, y_g, l_g)
      if cur_node in edge_node_bot:
        x_support = np.linspace(x_g - l_g, x_g + l_g, 100)
        y_support = np.linspace(y_g, y_g + l_g, 50)
        vertices = compute_S(x_support, y_support, x_g, y_g, l_g)
      if cur_node in edge_node_top:
        x_support = np.linspace(x_g - l_g, x_g + l_g, 100)
        y_support = np.linspace(y_g - l_g, y_g, 50)
        vertices = compute_S(x_support, y_support, x_g, y_g, l_g)
      if cur_node in edge_node_left:
        x_support = np.linspace(x_g, x_g + l_g, 50)
        y_support = np.linspace(y_g - l_g, y_g + l_g, 100)
        vertices = compute_S(x_support, y_support, x_g, y_g, l_g)
      if cur_node in edge_node_right:
        x_support = np.linspace(x_g - l_g, x_g, 50)
        y_support = np.linspace(y_g - l_g, y_g + l_g, 100)
        vertices = compute_S(x_support, y_support, x_g, y_g, l_g)
      if node not in boundary_nodes:
        x_support = np.linspace(x_g - l_g, x_g + l_g, 100)
        y_support = np.linspace(y_g - l_g, y_g + l_g, 100)
        vertices = compute_S(x_support, y_support, x_g, y_g, l_g)

  vertex_array = np.array(vertices)
  data = pv.PolyData(vertex_array)
  surf = data.delaunay_2d();
  return surf

def mpm_basis_linear(nx, node_coords, node):
  spacing = (node_coords[1][0] - node_coords[0][0])*0.5
  corner_node_mm = [0]
  corner_node_mp = [nx - 1]
  corner_node_pp = [nx*nx - 1]
  corner_node_pm = [nx*nx - nx]
  edge_node_bot = list(range(1, nx - 1))
  edge_node_top = list(range(nx*nx - nx + 1, nx*nx - 1))
  edge_node_left = list(range(nx, nx*nx - 2*nx + 1, nx))
  edge_node_right = list(range(2*nx-1, nx*nx - nx, nx))
  #edge_node_right = [7, 11]
  boundary_nodes = []
  boundary_nodes.extend(corner_node_mm)
  boundary_nodes.extend(corner_node_mp)
  boundary_nodes.extend(corner_node_pp)
  boundary_nodes.extend(corner_node_pm)
  boundary_nodes.extend(edge_node_bot)
  boundary_nodes.extend(edge_node_top)
  boundary_nodes.extend(edge_node_left)
  boundary_nodes.extend(edge_node_right)
  vertices = []
  for cur_node, point in enumerate(node_coords):
    if cur_node == node:
      vertices.append([point[0], point[1], 1.0])
      if cur_node in corner_node_mm:
        vertices.append([point[0]+spacing, point[1], 0.5])
        vertices.append([point[0]+spacing, point[1]+spacing, 0.25])
        vertices.append([point[0], point[1]+spacing, 0.5])
      if cur_node in corner_node_mp:
        vertices.append([point[0]-spacing, point[1]+spacing, 0.25])
        vertices.append([point[0]-spacing, point[1], 0.5])
        vertices.append([point[0], point[1]+spacing, 0.5])
      if cur_node in corner_node_pp:
        vertices.append([point[0]-spacing, point[1]-spacing, 0.25])
        vertices.append([point[0]-spacing, point[1], 0.5])
        vertices.append([point[0], point[1]-spacing, 0.5])
      if cur_node in corner_node_pm:
        vertices.append([point[0]+spacing, point[1]-spacing, 0.25])
        vertices.append([point[0]+spacing, point[1], 0.5])
        vertices.append([point[0], point[1]-spacing, 0.5])
      if cur_node in edge_node_bot:
        vertices.append([point[0]+spacing, point[1]+spacing, 0.25])
        vertices.append([point[0]+spacing, point[1], 0.5])
        vertices.append([point[0], point[1]+spacing, 0.5])
        vertices.append([point[0]-spacing, point[1]+spacing, 0.25])
        vertices.append([point[0]-spacing, point[1], 0.5])
      if cur_node in edge_node_top:
        vertices.append([point[0]+spacing, point[1]-spacing, 0.25])
        vertices.append([point[0]+spacing, point[1], 0.5])
        vertices.append([point[0], point[1]-spacing, 0.5])
        vertices.append([point[0]-spacing, point[1]-spacing, 0.25])
        vertices.append([point[0]-spacing, point[1], 0.5])
      if cur_node in edge_node_left:
        vertices.append([point[0]+spacing, point[1]-spacing, 0.25])
        vertices.append([point[0]+spacing, point[1]+spacing, 0.25])
        vertices.append([point[0], point[1]-spacing, 0.5])
        vertices.append([point[0], point[1]+spacing, 0.5])
        vertices.append([point[0]+spacing, point[1], 0.5])
      if cur_node in edge_node_right:
        vertices.append([point[0]-spacing, point[1]-spacing, 0.25])
        vertices.append([point[0]-spacing, point[1]+spacing, 0.25])
        vertices.append([point[0], point[1]-spacing, 0.5])
        vertices.append([point[0], point[1]+spacing, 0.5])
        vertices.append([point[0]-spacing, point[1], 0.5])
      if cur_node not in boundary_nodes:
        vertices.append([point[0]-spacing, point[1], 0.5])
        vertices.append([point[0]-spacing, point[1]-spacing, 0.25])
        vertices.append([point[0], point[1]-spacing, 0.5])
        vertices.append([point[0]+spacing, point[1]-spacing, 0.25])
        vertices.append([point[0]+spacing, point[1], 0.5])
        vertices.append([point[0]+spacing, point[1]+spacing, 0.25])
        vertices.append([point[0], point[1]+spacing, 0.5])
        vertices.append([point[0]-spacing, point[1]+spacing, 0.25])
    else:
      vertices.append([point[0], point[1], 0.0])
  vertex_array = np.array(vertices)

  data = pv.PolyData(vertex_array)
  surf = data.delaunay_2d();
  #print(surf)
  return surf
     
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
  
nx = 3
xnode, ynode, points = create_grid(0.0, 3.0, nx-1)
balls = create_node_balls(points)
particles, xpart, ypart = create_particle_balls(0.2, 2.8, nx*nx)
#shape_function()

#tab_cmap = plt.cm.get_cmap("RdYlGn", nx*nx)
tab_cmap = plt.cm.get_cmap("autumn_r", nx*nx)

pv.set_plot_theme('document')
p = pv.Plotter(notebook=False)

opacity = 1.0;
color = 0;
for i, point in enumerate(points):
  color = color + 1.0/(nx*nx)
  #opacity = opacity - 0.05
  #if (opacity < 0.3):
  # opacity = 1
  if (color >= 1.0):
    color = 0
  #grid = mpm_basis_linear(nx, points, i)
  grid = mpm_basis_actual(nx, nx, points, i)
  p.add_mesh(grid, color=tab_cmap(color)[0:3], opacity=opacity, specular=1.0, smooth_shading=False, show_edges=False)
  
for ball in balls:
  p.add_mesh(ball, color="red", show_edges=False)

for particle in particles:
  p.add_mesh(particle, color="green", show_edges=False)

p.add_axes()
p.add_bounding_box()
p.show_grid()
p.show()

