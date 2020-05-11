import numpy as np
import pyvista as pv

def mohr_coulomb(c, phi, low, high, numcells):

    rad = np.linspace(0.01, high, numcells)
    grid = pv.CylinderStructured(radius = rad, height = np.sqrt(3)*(high - low),
                                 center = [0.5*(high+low), 0.5*(high-low), 0.5*(high-low)],
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
  
c = 1.0e4
phi = 30

grid = mohr_coulomb(c, phi, -1.0e4, 1.0e6, 100)
contours = grid.contour([0.0])

df_dsigma = np.array([0.316228, 0, -0.948683])
stress_init = np.array([294450, 1.28206e-12, 1.2255e-12, 1.28206e-12, 86602.9, -5.02731e-32, 1.2255e-12, -5.02731e-32, 86602.9])
stress_trial = np.array([299422, 1.31845e-12, 1.31741e-12, 1.31845e-12, 88065.2, -5.16843e-32, 1.31741e-12, -5.16843e-32, 88065.2])
stress_init_eig, vec = np.linalg.eigh(np.reshape(stress_init, (3,3)))
stress_trial_eig, vec = np.linalg.eigh(np.reshape(stress_trial, (3,3)))
stress_init_eig = np.flip(stress_init_eig)
stress_trial_eig = np.flip(stress_trial_eig)

print(stress_init_eig)
stress_init_ball = pv.Sphere(radius=2000, center = stress_init_eig)
stress_trial_ball = pv.Sphere(radius=2000, center = stress_trial_eig)
df_dsigma_cyl, df_dsigma_cone = drawArrow(stress_init_eig, df_dsigma, 1.0e5, 1000)

stress_dir = stress_trial_eig - stress_init_eig
sigma_cyl, sigma_cone = drawArrow(stress_init_eig, stress_dir/np.linalg.norm(stress_dir), 1.0e5, 1000)

pv.set_plot_theme('document')
p = pv.Plotter(notebook=False)
p.add_mesh(grid, opacity=0.1)
p.add_mesh(contours, opacity=0.5, show_scalar_bar=False)
p.add_mesh(stress_init_ball, color="red", show_edges=False)
p.add_mesh(stress_trial_ball, color="red", show_edges=False)
p.add_mesh(df_dsigma_cyl, color="green", show_edges=False)
p.add_mesh(df_dsigma_cone, color="green", show_edges=False)
p.add_mesh(sigma_cyl, color="blue", show_edges=False)
p.add_mesh(sigma_cone, color="blue", show_edges=False)
p.add_axes()
p.add_bounding_box()
p.show_grid()
p.show()
