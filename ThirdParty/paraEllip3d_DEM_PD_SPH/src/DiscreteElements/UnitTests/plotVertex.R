# Requires:
#  install.packages('cda')
#  install.packages('rgl')
#  install.packages('orientlib')

require('cda')
require('rgl')
require('orientlib')

rm(list = ls())

#------------------------------------------------------------------
# Check if a matrix is a valid rotation matrix.
#------------------------------------------------------------------
isRotationMatrix <- function(R) {

  Rt <- t(R)
  shouldBeIdentity <- Rt %*% R
  I <- diag(3)
  n <- norm(I - shouldBeIdentity, "2")
  #print(n)
  return(n < 5e-6)
}
 
#------------------------------------------------------------------
# Calculates a x b 
#------------------------------------------------------------------
vec_cross <- function(aa, bb) {
  i1 <- c(2,3,1)
  i2 <- c(3,1,2)
  return (aa[i1]*bb[i2] - aa[i2]*bb[i1])
}

# Calculates rotation matrix to euler angles
# R[i,j] = E_i . e_j
#------------------------------------------------------------------
rotationMatrixToEulerAngles <- function(R) {
 
  #R12 = R[1,] %*% R[2,]
  #R13 = R[1,] %*% R[3,]
  #R23 = R[2,] %*% R[3,]

  #v3 = vec_cross(R[1,], R[2,])
  #v1 = vec_cross(R[2,], R[3,])
  #R[3,] = v3
  #v2 = vec_cross(R[3,], R[1,])
  #R[2,] = v2

  #R12 = R[1,] %*% R[2,]
  #R13 = R[1,] %*% R[3,]
  #R23 = R[2,] %*% R[3,]
  #print(paste("R12=", R12, "R13=", R13, "R23=", R23))

  #print("R:")
  #print(R)
  costheta = R[3,3]
  sinthetaSq = 1 - costheta^2
  if (sinthetaSq < 1.0e-20) {
    theta = acos(costheta)
    phi = atan2(R[2,1], R[1,1])
    psi = 0
  } else {
    theta = acos(costheta)
    phi = atan2(R[1,3], -R[2,3])
    psi = atan2(R[3,1],  R[3,2])
  }

  Rmat <- rotmatrix(t(R))
  angles <- eulerzxz(Rmat)
  phi_a = angles@x[1,1]
  theta_a = angles@x[1,2]
  psi_a = angles@x[1,3]

  #print(paste("phi:", phi, "theta:", theta, "psi:", psi))
  #print(paste("phi:", phi_a, "theta:", theta_a, "psi:", psi_a))
  #print(paste("phi:", phi_a - phi, "theta:", theta_a - theta, "psi:", psi_a - psi))

  cp = cos(phi)
  sp = sin(phi)
  ct = cos(theta)
  st = sin(theta)
  cs = cos(psi)
  ss = sin(psi)

  D = matrix(c(cp, -sp, 0, sp, cp, 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  C = matrix(c(1, 0, 0, 0 , ct, -st, 0, st, ct), nrow = 3, ncol = 3, byrow = TRUE)
  B = matrix(c(cs, -ss, 0, ss, cs, 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  DCB = D %*% C %*% B
  if (!isRotationMatrix(R) || !isRotationMatrix(DCB)) {
    diff = DCB - R 
    print(diff)
  }
  return(c(phi, theta, psi))
}

plotCDAEllipsoid <- function(index, pos, size, euler, dircos) {
  va = dircos[1:3]
  vb = dircos[4:6]
  vc = dircos[7:9]
  ppa = pos + size[1]*va
  ppb = pos + size[2]*vb
  ppc = pos + size[3]*vc
  arrow3d(pos, ppa, type="flat", col="red")
  arrow3d(pos, ppb, type="flat", col="green")
  arrow3d(pos, ppc, type="flat", col="blue")

  ell1 = rgl.ellipsoid(x = pos[1], y = pos[2], z = pos[3],
                       a = size[1], b = size[2], c = size[3],
                       phi = euler[1], theta = euler[2], psi = euler[3],
                       subdivide = 3, smooth = TRUE, color = "gold", alpha=0.5)
  plot3d(ell1, add = TRUE)
  text3d(x = pos[1], y = pos[2], z = pos[3], as.character(index))
}

my.ellipsoid <- function (pos, size, euler, subdivide = 3, smooth = TRUE, ...)
{
  x = pos[1]
  y = pos[2]
  z = pos[3]
  a = size[1]
  b = size[2]
  c = size[3]
  phi = euler[1]
  theta = euler[2]
  psi = euler[3]
  cp = cos(phi)
  sp = sin(phi)
  ct = cos(theta)
  st = sin(theta)
  cs = cos(psi)
  ss = sin(psi)
  D = matrix(c(cp, -sp, 0, sp, cp, 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  C = matrix(c(1, 0, 0, 0 , ct, -st, 0, st, ct), nrow = 3, ncol = 3, byrow = TRUE)
  B = matrix(c(cs, -ss, 0, ss, cs, 0, 0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  rotMat = D %*% C %*% B

  sphere <- rgl::subdivision3d(rgl::cube3d(...), subdivide)
  class(sphere) <- c("mesh3d","shape3d")

  norm <- sqrt(sphere$vb[1, ]^2 + 
               sphere$vb[2, ]^2 +
               sphere$vb[3, ]^2 )
  for (i in 1:3) sphere$vb[i, ] <- sphere$vb[i, ]/norm
  sphere$vb[4, ] <- 1
  sphere$normals <- sphere$vb
  result <- rgl::scale3d(sphere, a,b,c)
  rotM <- rotMat
  #rotM <- cpp_euler_passive(phi,theta,psi)
  result <- rgl::rotate3d(result,matrix=rotM)
  result <- rgl::translate3d(result, x,y,z)
  invisible(result)
}

my.ellipsoids <- function(positions, sizes, angles, colour = "red", ...){

  N <- NCOL(positions)
  colours <- rep(colour, length.out=N)
  ll <- lapply(seq(1,N), function(ii)
    my.ellipsoid(positions[,ii], sizes[,ii], angles[,ii],, col = colours[ii], ...))

  a1 <- arrow3d(c(0.002, 0, 0), c(0.0025, 0, 0), s = 1/3, width = 1/3, n = 10, type="rotation", col = "red", plot = FALSE)
  a2 <- arrow3d(c(0.002, 0, 0), c(0.002, 0.0005, 0), s = 1/3, width = 1/3, n = 10, type="rotation", plot = FALSE, col="green")
  a3 <- arrow3d(c(0.002, 0, 0), c(0.002, 0, 0.0005), s = 1/3, width = 1/3, n = 10, type="rotation", plot = FALSE, col="blue")
  rgl::shapelist3d(c(ll, list(a1, a2, a3)))

}

plotMyEllipsoid <- function(index, pos, size, euler, dircos) {
  va = dircos[1:3]
  vb = dircos[4:6]
  vc = dircos[7:9]
  ppa = pos + size[1]*va
  ppb = pos + size[2]*vb
  ppc = pos + size[3]*vc
  arrow3d(pos, ppa, type="flat", col="red")
  arrow3d(pos, ppb, type="flat", col="green")
  arrow3d(pos, ppc, type="flat", col="blue")

  ell2 = my.ellipsoid(pos, size, euler, subdivide = 3, smooth = TRUE, color = "gold", alpha=0.5)
  plot3d(ell2, add = TRUE)
  text3d(x = pos[1], y = pos[2], z = pos[3], as.character(index))
}

#------------------------------------------------------------------
# 1) Read the particle data
#------------------------------------------------------------------
#partDistCSV <- "input_particle_file_orig"
#partDistCSV <- "input_particle_file"
partDistCSV <- "particle_00001.dat"

# Read the first line (contains the number of particles)
df_numPart <- read.csv(partDistCSV, header = FALSE, nrows = 1, sep="",
                   blank.lines.skip = TRUE)
numPart <- df_numPart$V1

# Read the particle data into a data frame
# Data should contain the following header:
# > names(df)
# [1] "id"         "type"       "radius_a"   "radius_b"   "radius_c"  
# [6] "position_x" "position_y" "position_z" "axle_a_x"   "axle_a_y"  
#[11] "axle_a_z"   "axle_b_x"   "axle_b_y"   "axle_b_z"   "axle_c_x"  
#[16] "axle_c_y"   "axle_c_z"   "velocity_x" "velocity_y" "velocity_z"
#[21] "omga_x"     "omga_y"     "omga_z"     "force_x"    "force_y"   
#[26] "force_z"    "moment_x"   "moment_y"   "moment_z"  
df_part <- read.csv(partDistCSV, header = TRUE, skip = 1, sep="",
                   blank.lines.skip = TRUE, nrows = numPart)

#------------------------------------------------------------------
# 2) Create positions
#------------------------------------------------------------------
positions <- data.frame(x = df_part$position_x, y = df_part$position_y, z = df_part$position_z)
pos_mat <- data.matrix(t(positions))

#------------------------------------------------------------------
# 2) Create radii
#------------------------------------------------------------------
sizes <- data.frame(x = df_part$radius_a, y = df_part$radius_b, z = df_part$radius_c)
size_mat <- data.matrix(t(sizes))

#------------------------------------------------------------------
# 3) Create Euler angles
#------------------------------------------------------------------
angles <- data.frame(xa = df_part$axle_a_x, ya = df_part$axle_a_y, za = df_part$axle_a_z,
                     xb = df_part$axle_b_x, yb = df_part$axle_b_y, zb = df_part$axle_b_z,
                     xc = df_part$axle_c_x, yc = df_part$axle_c_y, zc = df_part$axle_c_z)
cos_angles <- cos(angles)
cos_mat <- data.matrix(cos_angles)
euler_angles <- apply(cos_angles, 1, 
                 function(x) {
                   #print(x)
                   R <- matrix(x, nrow = 3, ncol=3, byrow=TRUE)
                   #print(R)
                   angles = rotationMatrixToEulerAngles(R)
                   return(angles)
                 })
euler_mat <- data.matrix(euler_angles)

#------------------------------------------------------------------
# 4) Plot arrows
#------------------------------------------------------------------
#open3d()
#arrow3d(c(0, 0, 0), c(1, 0, 0), s = 1/7, width = 1/5, n = 10, type="rotation", col="red")
#arrow3d(c(0, 0, 0), c(0, 1, 0), s = 1/7, width = 1/5, n = 10, type="rotation", col="green")
#arrow3d(c(0, 0, 0), c(0, 0, 1), s = 1/7, width = 1/5, n = 10, type="rotation", col="blue")
#axes3d()

#------------------------------------------------------------------
# 5) Plot rgl ellipsoids
#------------------------------------------------------------------
#open3d()
#rgl.ellipsoids(pos_mat, size_mat, euler_angles, col="gold", add=TRUE)
#axes3d()
#view_matrix = matrix(c(-0.21751294, -0.53017938,  0.81950366,  0.00000000,  
#                        0.97544736, -0.08853958,  0.20162290,  0.00000000,
#                       -0.03433822,  0.84324223,  0.53642339,  0.00000000,
#                        0.00000000,  0.00000000,  0.00000000,  1.00000000),
#                     nrow = 4, ncol = 4)
#view3d(userMatrix = view_matrix)
#rgl.snapshot("rgl_cda_ellipsoids.png")

#------------------------------------------------------------------
# 6) Plot individual rgl ellipsoids
#------------------------------------------------------------------
#open3d()
#for (i in c(2,3)) {
  #print(euler_mat[,i])
#  plotCDAEllipsoid(i, pos_mat[,i], size_mat[,i], euler_mat[,i], cos_angles[i,])
#}
#arrow3d(c(0, 0, 0), c(0.0006, 0, 0), s = 1/7, width = 1/5, n = 10, type="rotation", col="red")
#arrow3d(c(0, 0, 0), c(0, 0.0006, 0), s = 1/7, width = 1/5, n = 10, type="rotation", col="green")
#arrow3d(c(0, 0, 0), c(0, 0, 0.0006), s = 1/7, width = 1/5, n = 10, type="rotation", col="blue")
#axes3d()
#view3d(fov=0)

#------------------------------------------------------------------
# 7) Plot individual my ellipsoids
#------------------------------------------------------------------
#open3d()
#for (i in c(2,3)) {
#  #print(euler_mat[,i])
#  plotMyEllipsoid(i, pos_mat[,i], size_mat[,i], euler_mat[,i], cos_angles[i,])
#}
#arrow3d(c(0, 0, 0), c(0.0006, 0, 0), s = 1/7, width = 1/5, n = 10, type="rotation", col="red")
#arrow3d(c(0, 0, 0), c(0, 0.0006, 0), s = 1/7, width = 1/5, n = 10, type="rotation", col="green")
#arrow3d(c(0, 0, 0), c(0, 0, 0.0006), s = 1/7, width = 1/5, n = 10, type="rotation", col="blue")
#axes3d()
#view3d(fov=0)

#------------------------------------------------------------------
# Plot the bounding boxes
#------------------------------------------------------------------
#open3d()
#cube <- cube3d(color = "green", alpha = 0.3)

#vertices = matrix(c(-10,-10,-10, 1,  10,-10,-10, 1, 
#                     10, 12,-10, 1, -10, 12,-10, 1,
#                    -10,-10, 14, 1,  10,-10, 14, 1, 
#                     10, 12, 14, 1, -10, 12, 14, 1), 
#                  nrow=8, byrow = TRUE)
#xmax = 10
#ymax = 10
#zmax = 10
#xmax = 10
#ymax = 12
#zmax = 14
#cube_scaled <- rgl::scale3d(cube, xmax, ymax, zmax)
#shade3d(cube_scaled)
#points3d(t(cube_scaled$vb))
#texts3d(vertices, adj = c(1.5, 1.5), texts=c(1,2,3,4,5,6,7,8))
#for (i in 1:6)
#    lines3d(t(cube_scaled$vb)[cube_scaled$ib[,i],])
#
#arrow3d(c(0, 0, 0), c(xmax, 0, 0), s = 1/7, width = 1/5, n = 10, type="rotation", col="red")
#arrow3d(c(0, 0, 0), c(0, ymax, 0), s = 1/7, width = 1/5, n = 10, type="rotation", col="green")
#arrow3d(c(0, 0, 0), c(0, 0, zmax), s = 1/7, width = 1/5, n = 10, type="rotation", col="blue")
#texts3d(c(1.05*xmax, 0, 0), adj = 1, texts=c("x"))
#texts3d(c(0, 1.05*ymax, 0), adj = 1, texts=c("y"))
#texts3d(c(0, 0, 1.05*zmax), adj = 1, texts=c("z"))

#------------------------------------------------------------------
# Plot the faces with different colors
#------------------------------------------------------------------
#zoom<-par3d()$zoom
#userMatrix<-par3d()$userMatrix
#windowRect<-par3d()$windowRect

zoom = 1.0
userMatrix = matrix(c( 0.6981224, -0.7159739, -0.0007422371, 0,
                      -0.3004532, -0.2939027,  0.9073821902, 0,
                      -0.6498820, -0.6332429, -0.4202980101, 0,
                       0.0000000,  0.0000000,  0.0000000000, 1), 4, 4,
                       byrow = TRUE)
windowRect = c(152, 164, 1103, 966)

open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
cube1 <- cube3d()
vertices = matrix(c(-10,-10,-10, 1,  10,-10,-10, 1, 
                     10, 10,-10, 1, -10, 10,-10, 1,
                    -10,-10, 10, 1,  10,-10, 10, 1, 
                     10, 10, 10, 1, -10, 10, 10, 1), 
                  nrow=8, byrow = TRUE)
xmax = 10
ymax = 10
zmax = 10
cube1_scaled <- rgl::scale3d(cube1, xmax, ymax, zmax)
shade3d(cube1_scaled, 
        col = c(rep(4,4), rep(8, 4), rep(8,4), rep(2,4), rep(3,4), rep(8,4)),
        alpha = 0.7)
points3d(t(cube1_scaled$vb))
#texts3d(vertices, adj = c(2, 2), texts=c(1,2,3,4,5,6,7,8))
for (i in 1:6)
    lines3d(t(cube1_scaled$vb)[cube1_scaled$ib[,i],])

#------------------------------------------------------------------
# Plot the face normal arrows
#------------------------------------------------------------------
#arrow3d(c(-xmax, 0, 0), c(-2*xmax, 0, 0), barblen = 0.1, s = 0.5, width = 0.3, n = 10, type="rotation", col="#ff5733")
#arrow3d(c(0, -ymax, 0), c(0, -ymax-xmax, 0), barblen = 0.1, s = 0.5, width = 0.3, n = 10, type="rotation", col="#54f576")
#arrow3d(c(0, 0, -zmax), c(0, 0, -zmax-xmax), barblen = 0.1, s = 0.5, width = 0.3, n = 10, type="rotation", col="#54bdf5")
#texts3d(c(-2*xmax, 0, 0), adj = c(1.5,1.5), texts=c("x-"))
#texts3d(c(0, -ymax-xmax, 0), adj = c(1.5,1.5), texts=c("y-"))
#texts3d(c(0, 0, -zmax-xmax), adj = c(1.5,1.5), texts=c("z-"))

#------------------------------------------------------------------
# Plot the circular arrows
#------------------------------------------------------------------
# Parametric curves
#theta <- seq(from = 270, to = 0, length = 10)
#rr <- 5
#xs <- rr*cos(theta*pi/180)
#ys <- rr*sin(theta*pi/180)
#zs <- rr*sin(theta*0) - 1.2*zmax
#tt <- seq(from = 0, to = 3, length = 10)
#df1 = data.frame(t = tt, x = xs, y = ys, z = zs)
#ts <- seq( from = min(df1$t), max(df1$t), length=100 )
#
#df2 = apply( df1[,-1], 2, function(u) spline( df1$t, u, xout = ts )$y )
#lines3d( df2, col="yellow", lwd=5)
#vv = df2[100,] - df2[90,]
#arrow3d(df2[95,], df2[95,]+vv, barblen = 0.05, s = 0.5, width = 0.3, n = 10, type="rotation", col="yellow")
#
#df3 = data.frame(x = df2[,1], y = df2[,3], z = df2[,2])
#lines3d( df3, col="yellow", lwd=5)
#vv = df3[1,] - df3[10,]
#arrow3d(df3[5,], df3[5,]+vv, barblen = 0.05, s = 0.5, width = 0.3, n = 10, type="rotation", col="yellow")
#
#df4 = data.frame(x = df2[,3], y = df2[,2], z = df2[,1])
#lines3d( df4, col="yellow", lwd=5)
#vv = df4[1,] - df4[10,]
#arrow3d(df4[5,], df4[5,]+vv, barblen = 0.05, s = 0.5, width = 0.3, n = 10, type="rotation", col="yellow")


#------------------------------------------------------------------
# Plot the extended cube
#------------------------------------------------------------------
max_size = 3.6
cube_extra <- cube3d(color = "blue", alpha = 0.1)
cube_extra_scaled <- rgl::scale3d(cube_extra, xmax + max_size, ymax + max_size, zmax + max_size)
cube_extra_shift <- rgl::translate3d(cube_extra_scaled, 1.01*max_size, 1.01*max_size, 1.01*max_size)
shade3d(cube_extra_shift)
for (i in 1:6)
    lines3d(t(cube_extra_shift$vb)[cube_extra_shift$ib[,i],])


#------------------------------------------------------------------
# Plot my ellipsoids (vertex 1)
#------------------------------------------------------------------
#ids <- df_part$id
#my.ellipsoids(pos_mat[,c(2,9)], size_mat[,c(2,9)], euler_angles[,c(2,9)], col="gold", add=TRUE, alpha = 1.0)
#texts3d(t(pos_mat)[c(2,9),], adj = c(1.5, 1.5), texts=ids[c(2,9)])
#arrow3d(pos_mat[,2], pos_mat[,9], s = 0.07, width = 0.1, n = 10, type="rotation", col="blue")

#------------------------------------------------------------------
# Plot my ellipsoids (vertex 2)
#------------------------------------------------------------------
#ids <- df_part$id
#my.ellipsoids(pos_mat[,c(5,14,15,16)], size_mat[,c(5,14,15,16)], euler_angles[,c(5,14,15,16)], col="gold", add=TRUE, alpha = 1.0)
#texts3d(t(pos_mat)[c(5,14,15,16),], adj = c(1.5, 1.5), texts=ids[c(5,14,15,16)])
#arrow3d(pos_mat[,5], pos_mat[,14], s = 0.07, width = 0.1, n = 10, type="rotation", col="blue")
#arrow3d(pos_mat[,5], pos_mat[,15], s = 0.07, width = 0.1, n = 10, type="rotation", col="blue")
#arrow3d(pos_mat[,5], pos_mat[,16], s = 0.07, width = 0.1, n = 10, type="rotation", col="blue")
#rgl.snapshot("rgl_my_ellipsoids.png")

#------------------------------------------------------------------
# Plot my ellipsoids (vertex 3)
#------------------------------------------------------------------
ids <- df_part$id
my.ellipsoids(pos_mat[,c(8,21:27)], size_mat[,c(8,21:27)], euler_angles[,c(8,21:27)], col="gold", add=TRUE, alpha = 1.0)
texts3d(t(pos_mat)[c(8,21:27),], adj = c(1.5, 1.5), texts=ids[c(8,21:27)])
arrow3d(pos_mat[,8], pos_mat[,21], s = 0.07, width = 0.1, n = 10, type="rotation", col="blue")
arrow3d(pos_mat[,8], pos_mat[,22], s = 0.07, width = 0.1, n = 10, type="rotation", col="blue")
arrow3d(pos_mat[,8], pos_mat[,23], s = 0.07, width = 0.1, n = 10, type="rotation", col="blue")
arrow3d(pos_mat[,8], pos_mat[,24], s = 0.07, width = 0.1, n = 10, type="rotation", col="blue")
arrow3d(pos_mat[,8], pos_mat[,25], s = 0.07, width = 0.1, n = 10, type="rotation", col="blue")
arrow3d(pos_mat[,8], pos_mat[,26], s = 0.07, width = 0.1, n = 10, type="rotation", col="blue")
arrow3d(pos_mat[,8], pos_mat[,27], s = 0.07, width = 0.1, n = 10, type="rotation", col="blue")
#rgl.snapshot("rgl_my_ellipsoids.png")

#ids <- df_part$id
#my.ellipsoids(pos_mat[,c(-1, -3, -4, -7, -10,-11,-12,-13,-6,-17,-18,-19,-20,-2,-5,-9,-14,-15,-16)], size_mat[,c(-1, -3, -4, -7, -10,-11,-12,-13,-6,-17,-18,-19,-20,-2,-5,-9,-14,-15,-16)], euler_angles[,c(-1, -3, -4, -7, -10,-11,-12,-13,-6,-17,-18,-19,-20,-2,-5,-9,-14,-15,-16)], col="gold", add=TRUE, alpha = 1.0)
#texts3d(t(pos_mat)[c(-1, -3, -4, -7, -10,-11,-12,-13,-6,-17,-18,-19,-20,-2,-5,-9,-14,-15,-16),], adj = c(1.5, 1.5), texts=ids[c(-1, -3, -4, -7, -10,-11,-12,-13,-6,-17,-18,-19,-20,-2,-5,-9,-14,-15,-16)])

#ids <- df_part$id
#my.ellipsoids(pos_mat, size_mat, euler_angles, col="gold", add=TRUE, alpha = 1.0)
#texts3d(t(pos_mat), adj = c(1.5, 1.5), texts=ids)

