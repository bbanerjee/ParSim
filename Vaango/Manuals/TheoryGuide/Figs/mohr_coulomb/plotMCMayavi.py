import numpy as np
from skimage import measure
from mayavi import mlab

def mohr_coulomb(c, phi, low, high, numcells):

    sig = np.linspace(low, high, numcells)

    arr = np.array([[[1.0 for x in range(sig.shape[0])]
                     for y in range(sig.shape[0])]
                    for z in range(sig.shape[0])])

    sin_phi = np.sin(phi * np.pi / 180.0)
    cos_phi = np.cos(phi * np.pi / 180.0)

    for i, s1 in enumerate(sig):
        for j, s2 in enumerate(sig):
            for k, s3 in enumerate(sig):
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

    return arr

def plot_mayavi(verts_mc, faces_mc):
    mlab.triangular_mesh(
        [vert[0] for vert in verts_mc],
        [vert[1] for vert in verts_mc],
        [vert[2] for vert in verts_mc],
        faces_mc,
    )

    mlab.show()


c = 0.5
phi = 30

surf_arr = mohr_coulomb(c, phi, -1, 3, 100)

verts_mc, faces_mc, normals, values = measure.marching_cubes_lewiner(
    surf_arr, 0)

plot_mayavi(verts_mc, faces_mc)


#import matplotlib.pyplot as plt
#import visvis as vv
#from skimage.draw import ellipsoid
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection

#def plot_vv():
#    vv.mesh(np.fliplr(verts_mc), faces_mc, normals, values)
#    vv.use().Run()


#def plot_mp():
#    mesh_mc = Poly3DCollection(verts_mc[faces_mc])
#    # mesh_mc.set_edgecolor('k')
#
#    fig = plt.figure(figsize=(10, 10))
#    ax = fig.add_subplot(111, projection="3d")
#    ax.add_collection3d(mesh_mc)
#    ax.set_xlabel("sigma_1")
#    ax.set_ylabel("sigma_2")
#    ax.set_zlabel("sigma_3")
#
#    # ax.set_xlim(-1, 3)
#    # ax.set_ylim(-1, 3)
#    # ax.set_zlim(-1, 3)
#
#    ax.set_xlim(0, 30)
#    ax.set_ylim(0, 30)
#    ax.set_zlim(0, 30)
#
#    plt.tight_layout()
#    plt.show()
