import numpy as np
import pyvista as pv

def rotationMatrix(axis, angle):
  axis = axis / np.linalg.norm(axis)
  I = np.identity(3, dtype=float)
  A = np.array([0, -axis[2], axis[1], axis[2], 0, -axis[0], -axis[1], axis[0], 0]).reshape(3,3)
  axis = axis.reshape(1,3)
  aa = np.matmul(axis.transpose(), axis)
  R = (I - aa)*np.cos(angle) + aa - A*np.sin(angle)
  #print("angle = ", angle * 180 / np.pi)
  #print("axis = ", axis)
  #print("aa = ", aa)
  #print("A = ", A)
  #print("R = ", R)
  #print("RRT = ", np.matmul(R, R.transpose()))
  #print("RTR = ", np.matmul(R.transpose(), R))
  return R

def drawArrow(start, direction, scale):
  length = np.linalg.norm(direction)
  dir_hat = direction / length
  cyl_bot = start
  cyl_top = start + dir_hat * length
  cyl_cen = 0.5 * (cyl_bot + cyl_top)
  arrow_cyl = pv.Cylinder(cyl_cen, direction = dir_hat, radius = scale, 
                          height = length)
  cone_bot = cyl_top
  cone_top = cone_bot + direction * length * 0.3
  cone_cen = 0.5 * (cone_bot + cone_top)
  arrow_cone = pv.Cone(cone_cen, dir_hat, height = np.linalg.norm(cone_top - cone_bot), 
                       radius = 2.5*scale)

  return arrow_cyl, arrow_cone

def drawLine(start, direction, scale):
  length = np.linalg.norm(direction)
  dir_hat = direction / length
  cyl_bot = start
  cyl_top = start + dir_hat * length
  cyl_cen = 0.5 * (cyl_bot + cyl_top)
  cyl = pv.Cylinder(cyl_cen, direction = dir_hat, radius = scale, 
                    height = length)
  return cyl

def drawCube(center, v1, v2, v3, scale):
  start = center - v1 - v2 - v3 
  l1 = drawLine(start, 2.0 * v1, scale)
  l2 = drawLine(start, 2.0 * v2, scale)
  l3 = drawLine(start, 2.0 * v3, scale)
  start = center + v1 + v2 + v3 
  l4 = drawLine(start, -2.0 * v1, scale)
  l5 = drawLine(start, -2.0 * v2, scale)
  l6 = drawLine(start, -2.0 * v3, scale)
  start = center - v1 - v2 + v3 
  l7 = drawLine(start, 2.0 * v1, scale)
  l8 = drawLine(start, 2.0 * v2, scale)
  start = center - v1 + v2 + v3 
  l9 = drawLine(start, -2.0 * v3, scale)
  start = center + v1 + v2 - v3 
  l10 = drawLine(start, -2.0 * v1, scale)
  l11 = drawLine(start, -2.0 * v2, scale)
  start = center + v1 - v2 - v3 
  l12 = drawLine(start, 2.0 * v3, scale)
  return (l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12)
  

  
center = np.array([0.0, 0.0, 0.0])
v1 = np.array([0.353553, 0.353553, 0.4])
v2_ref = np.array([-0.707107, 0.707107, 1]) 
v2 = np.cross(v1/np.linalg.norm(v1), v2_ref/np.linalg.norm(v2_ref)) 
v2 = v2/np.linalg.norm(v2) 
v3 = np.cross(v1/np.linalg.norm(v1), v2/np.linalg.norm(v2)) 
v3 = v3/np.linalg.norm(v3) * 2.0
v1_cyl, v1_cone = drawArrow(center, v1, 0.01)
v2_cyl, v2_cone = drawArrow(center, v2, 0.01)
v3_cyl, v3_cone = drawArrow(center, v3, 0.01)
np.set_printoptions(precision=12)
print("v1 = ", v1)
print("v2 = ", v2)
print("||v1|| = ", np.linalg.norm(v1))
print("||v2|| = ", np.linalg.norm(v2))
print("||v3|| = ", np.linalg.norm(v3))

x_axis = np.array([1, 0, 0])
rot_axis = np.cross(v1/np.linalg.norm(v1), x_axis)
rot_angle = np.arccos(np.dot(v1/np.linalg.norm(v1), x_axis))
R_x = rotationMatrix(rot_axis, rot_angle)
rot_cyl, rot_cone = drawArrow(center, rot_axis, 0.01)

#v1_rot = np.matmul(R_x, v1)
#v1_rot_cyl, v1_rot_cone = drawArrow(center, v1_rot, 0.01)
v1_rotT = np.matmul(R_x.transpose(), v1)
v1_rotT_cyl, v1_rotT_cone = drawArrow(center, v1_rotT, 0.01)
v2_rotT = np.matmul(R_x.transpose(), v2)
v2_rotT_cyl, v2_rotT_cone = drawArrow(center, v2_rotT, 0.01)
v3_rotT = np.matmul(R_x.transpose(), v3)
v3_rotT_cyl, v3_rotT_cone = drawArrow(center, v3_rotT, 0.01)

y_axis = np.array([0, 1, 0])
rot_axis = np.cross(v2_rotT/np.linalg.norm(v2_rotT), y_axis)
rot_angle = np.arccos(np.dot(v2_rotT/np.linalg.norm(v2_rotT), y_axis))
R_y = rotationMatrix(rot_axis, rot_angle)
rot_cyl, rot_cone = drawArrow(center, rot_axis, 0.01)

v1_rotT_y = np.matmul(R_y.transpose(), v1_rotT)
v1_rotT_y_cyl, v1_rotT_y_cone = drawArrow(center, v1_rotT_y, 0.01)
v2_rotT_y = np.matmul(R_y.transpose(), v2_rotT)
v2_rotT_y_cyl, v2_rotT_y_cone = drawArrow(center, v2_rotT_y, 0.01)
v3_rotT_y = np.matmul(R_y.transpose(), v3_rotT)
v3_rotT_y_cyl, v3_rotT_y_cone = drawArrow(center, v3_rotT_y, 0.01)

# Total rotation matrix
R_xy = np.matmul(R_y.transpose(), R_x.transpose())
v1_rotT_xy = np.matmul(R_xy, v1)
v1_rotT_xy_cyl, v1_rotT_xy_cone = drawArrow(center, v1_rotT_xy, 0.01)
v2_rotT_xy = np.matmul(R_xy, v2)
v2_rotT_xy_cyl, v2_rotT_xy_cone = drawArrow(center, v2_rotT_xy, 0.01)
v3_rotT_xy = np.matmul(R_xy, v3)
v3_rotT_xy_cyl, v3_rotT_xy_cone = drawArrow(center, v3_rotT_xy, 0.01)

# Draw bounding box
line1 = drawLine(center - v1, 2.0 * v1, 0.005)
line2 = drawLine(center - v2, 2.0 * v2, 0.005)
line3 = drawLine(center - v3, 2.0 * v3, 0.005)
l0 = drawCube(center, v1, v2, v3, 0.005)

# Draw bounding box
line1 = drawLine(center - v1_rotT_xy, 2.0 * v1_rotT_xy, 0.005)
line2 = drawLine(center - v2_rotT_xy, 2.0 * v2_rotT_xy, 0.005)
line3 = drawLine(center - v3_rotT_xy, 2.0 * v3_rotT_xy, 0.005)
l1 = drawCube(center, v1_rotT_xy, v2_rotT_xy, v3_rotT_xy, 0.005)

pv.set_plot_theme('document')
p = pv.Plotter(notebook=False)
#p.add_mesh(v1_cyl, color="red", show_edges=False)
#p.add_mesh(v1_cone, color="red", show_edges=False)
#p.add_mesh(v2_cyl, color="green", show_edges=False)
#p.add_mesh(v2_cone, color="green", show_edges=False)
#p.add_mesh(v3_cyl, color="blue", show_edges=False)
#p.add_mesh(v3_cone, color="blue", show_edges=False)
#p.add_mesh(rot_cyl, color="magenta", show_edges=False)
#p.add_mesh(rot_cone, color="magenta", show_edges=False)
#p.add_mesh(v1_rot_cyl, color="yellow", show_edges=False)
#p.add_mesh(v1_rot_cone, color="yellow", show_edges=False)
#p.add_mesh(v1_rotT_cyl, color="cyan", show_edges=False)
#p.add_mesh(v1_rotT_cone, color="cyan", show_edges=False)
#p.add_mesh(v2_rotT_cyl, color="#FF7E00", show_edges=False)
#p.add_mesh(v2_rotT_cone, color="#FF7E00", show_edges=False)
#p.add_mesh(v3_rotT_cyl, color="#BFFF00", show_edges=False)
#p.add_mesh(v3_rotT_cone, color="#BFFF00", show_edges=False)
#p.add_mesh(v1_rotT_y_cyl, color="cyan", show_edges=False)
#p.add_mesh(v1_rotT_y_cone, color="cyan", show_edges=False)
#p.add_mesh(v2_rotT_y_cyl, color="#FF7E00", show_edges=False)
#p.add_mesh(v2_rotT_y_cone, color="#FF7E00", show_edges=False)
#p.add_mesh(v3_rotT_y_cyl, color="#BFFF00", show_edges=False)
#p.add_mesh(v3_rotT_y_cone, color="#BFFF00", show_edges=False)
p.add_mesh(v1_rotT_xy_cyl, color="red", show_edges=False)
p.add_mesh(v1_rotT_xy_cone, color="red", show_edges=False)
p.add_mesh(v2_rotT_xy_cyl, color="green", show_edges=False)
p.add_mesh(v2_rotT_xy_cone, color="green", show_edges=False)
p.add_mesh(v3_rotT_xy_cyl, color="blue", show_edges=False)
p.add_mesh(v3_rotT_xy_cone, color="blue", show_edges=False)
#p.add_mesh(line1, color="black", show_edges=False)
#p.add_mesh(line2, color="black", show_edges=False)
#p.add_mesh(line3, color="black", show_edges=False)
for line in l0:
  p.add_mesh(line, color="black", show_edges=False)
for line in l1:
  p.add_mesh(line, color="cyan", show_edges=False)

p.add_axes()
#p.add_bounding_box()
#p.show_grid()
p.show()
