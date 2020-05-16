import numpy as np
import pyvista as pv

def mohr_coulomb(c, phi, low, high, numcells):

    sig = np.linspace(low, high, numcells)

    arr = np.array([[[1.0 for x in range(sig.shape[0])]
                     for y in range(sig.shape[0])]
                    for z in range(sig.shape[0])])

    xcoord, ycoord, zcoord = np.meshgrid(sig, sig, sig)

    sin_phi = np.sin(phi * np.pi / 180.0)
    cos_phi = np.cos(phi * np.pi / 180.0)

    for i, s1 in enumerate(sig):
        for j, s2 in enumerate(sig):
            for k, s3 in enumerate(sig):
                xcoord[i, j, k] = s1;
                ycoord[i, j, k] = s2;
                zcoord[i, j, k] = s3;
                sig0 = np.array([s1, s2, s3])
                ind = np.argsort(sig0)
                i1 = ind[2]
                i3 = ind[0]
                arr[i, j, k] = \
                  ((sig0[i1] - sig0[i3]) -
                  (sig0[i1] + sig0[i3]) * sin_phi - 
                  2.0 * c * cos_phi) / \
                  (abs(sig0[i1]) + abs(sig0[i3]) + 2.0 * c)
                if s1 + s2 + s3 > (high - low):
                    arr[i, j, k] = 100

    return arr, xcoord, ycoord, zcoord

c = 1.0e4
phi = 30

surf_arr, x, y, z = mohr_coulomb(c, phi, -1.0e4, 1.0e6, 150)
print(x.shape, surf_arr.shape)

grid = pv.StructuredGrid(x, y, z)
grid["data"] = np.reshape(surf_arr, -1)

contours = grid.contour([0.0])

pv.set_plot_theme('document')
p = pv.Plotter(notebook=False)
p.add_mesh(grid, opacity=0.1)
p.add_mesh(contours, scalars=contours.points[:, 2], show_scalar_bar=False,
           opacity=0.5)
p.add_axes()
p.add_bounding_box()
p.show_grid()
p.show()
