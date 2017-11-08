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
partDistCSV <- "particle_00000.dat"

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
# 8) Plot my ellipsoids
#------------------------------------------------------------------
open3d()
my.ellipsoids(pos_mat, size_mat, euler_angles, col="gold", add=TRUE)
axes3d()
view_matrix = matrix(c(-0.21751294, -0.53017938,  0.81950366,  0.00000000,  
                        0.97544736, -0.08853958,  0.20162290,  0.00000000,
                       -0.03433822,  0.84324223,  0.53642339,  0.00000000,
                        0.00000000,  0.00000000,  0.00000000,  1.00000000),
                     nrow = 4, ncol = 4)
view3d(userMatrix = view_matrix)
#rgl.snapshot("rgl_my_ellipsoids.png")

#------------------------------------------------------------------
# 9) Plot the bounding boxes
#------------------------------------------------------------------
cube <- cube3d(color = "green", alpha = 0.3)
xmin = c(-10, -10, -10)
xmax = c(10, 10, 10)
cube_scaled <- rgl::scale3d(cube, 10, 10, 10)
shade3d(cube_scaled)

max_size = 3.6
cube_extra <- cube3d(color = "blue", alpha = 0.1)
cube_extra_scaled <- rgl::scale3d(cube_extra, 10 + max_size, 10 + max_size, 10 + max_size)
cube_extra_shift <- rgl::translate3d(cube_extra_scaled, max_size, max_size, max_size)
shade3d(cube_extra_shift)
for (i in 1:6)
    lines3d(t(cube_extra_shift$vb)[cube_extra_shift$ib[,i],])

