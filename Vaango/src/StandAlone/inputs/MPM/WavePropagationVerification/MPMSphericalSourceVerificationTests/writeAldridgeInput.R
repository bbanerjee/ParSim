# Set the working directory and jre location
setwd(".")

# Open file for writing
con_hires = file("input_param_file_hires.in", "w+")
con_lores = file("input_param_file_lores.in", "w+")

# Get the mass density
#  <density> 1700 </density>
rho = 1700

# Get the bulk and shear moduli
#  <K> 6.0e7 </K>
#  <G> 9.0e7</G>
K = 6.0e7
G = 9.0e7
M = K + 4.0/3.0*G

# Compute alpha, beta
alpha = sqrt(M/rho)
beta = sqrt(G/rho)

# Write alpha, beta, rho
cat(alpha, beta, rho, "\n", file = con_hires, sep=",")
cat(alpha, beta, rho, "\n", file = con_lores, sep=",")

# Get tmin, dt, tlen
#  <maxTime>            2.00e-2     </maxTime>
#  <initTime>           0.0         </initTime>
#  <outputInterval> 1.0e-7 </outputInterval>
tmin = 0.0
dt = 5.0e-7
tlen = 2.5e-3

# Write tmin, dt, tlen
cat(tmin, dt, tlen, "\n", file = con_hires, sep=",")
cat(tmin, dt, tlen, "\n", file = con_lores, sep=",")

# Get sphere coordinates
#  <origin> [0.0, 0.00, 0.60] </origin>
#  <radius> 0.02 </radius>
xs = 0.0
ys = 0.0
zs = 0.60
rs = 0.02

# Write xs, ys, zs
cat(xs, ys, zs, "\n", file = con_hires, sep=",")
cat(xs, ys, zs, "\n", file = con_lores, sep=",")

# Write rs
cat(rs, "\n", file = con_hires, sep=",")
cat(rs, "\n", file = con_lores, sep=",")

# Source amplitude
# 1 x value specified in input file containing the load curve
samp = 1

# Write samp
cat(samp, "\n", file = con_hires, sep=",")
cat(samp, "\n", file = con_lores, sep=",")

# Source type
stype = 4

# Write stype
cat(stype, "\n", file = con_hires, sep=",")
cat(stype, "\n", file = con_lores, sep=",")

# Number of receivers
# 6 disp, 6 vel, 6 acceleration, 6 pressure
nr = 24
rtype1 = 0
rtype2 = 1
rtype3 = 2
rtype4 = 4

# Write nr
cat(nr, "\n", file = con_hires, sep=",")
cat(nr, "\n", file = con_lores, sep=",")

# Receiver coordinates
#  rad = [-0.1, -0.152, -0.203, -0.254, -0.305, -0.356]
rad_acc_hires = c(0.105, 0.155, 0.205, 0.255, 0.305, 0.355)
xc_acc_hires = 0.0
yc_acc_hires = 0.005
zc_acc_hires = 0.605

rad_acc_lores = c(0.11, 0.15, 0.21, 0.25, 0.31, 0.35)
xc_acc_lores = 0.0
yc_acc_lores = 0.01
zc_acc_lores = 0.61

rad_data_hires = c(0.0975, 0.1475, 0.1975, 0.2525, 0.3025, 0.3525)
xc_data_hires = 0.0
yc_data_hires = 0.0025
zc_data_hires = 0.6025

rad_data_lores = c(0.095, 0.155, 0.205, 0.255, 0.305, 0.355)
xc_data_lores = 0.0
yc_data_lores = 0.005
zc_data_lores = 0.605

# Compute coordinates
theta0 = 0
phi0 = 0

xr_acc_hires = xc_acc_hires + rad_acc_hires*cos(theta0*pi/180.0)*cos(phi0*pi/180.0)
yr_acc_hires = yc_acc_hires + rad_acc_hires*sin(theta0*pi/180.0)*cos(phi0*pi/180.0)
zr_acc_hires = zc_acc_hires + rad_acc_hires*sin(phi0*pi/180.0)
pos = cbind(xr_acc_hires, yr_acc_hires, zr_acc_hires)

xr_acc_lores = xc_acc_lores + rad_acc_lores*cos(theta0*pi/180.0)*cos(phi0*pi/180.0)
yr_acc_lores = yc_acc_lores + rad_acc_lores*sin(theta0*pi/180.0)*cos(phi0*pi/180.0)
zr_acc_lores = zc_acc_lores + rad_acc_lores*sin(phi0*pi/180.0)

xr_data_hires = xc_data_hires + rad_data_hires*cos(theta0*pi/180.0)*cos(phi0*pi/180.0)
yr_data_hires = yc_data_hires + rad_data_hires*sin(theta0*pi/180.0)*cos(phi0*pi/180.0)
zr_data_hires = zc_data_hires + rad_data_hires*sin(phi0*pi/180.0)

xr_data_lores = xc_data_lores + rad_data_lores*cos(theta0*pi/180.0)*cos(phi0*pi/180.0)
yr_data_lores = yc_data_lores + rad_data_lores*sin(theta0*pi/180.0)*cos(phi0*pi/180.0)
zr_data_lores = zc_data_lores + rad_data_lores*sin(phi0*pi/180.0)

# Amplitude at receiver
ramp = 1.0

# Polar and azimuthal angles
theta = 90.0 # From +z axis
phi = 0.0    # From +x axis

# Write displacement receiver data
for (ii in 1:nrow(pos)) {
  cat(rtype1, xr_data_hires[ii], yr_data_hires[ii], zr_data_hires[ii], ramp, theta, phi, "\n", 
      file = con_hires, sep = ",")
  cat(rtype1, xr_data_lores[ii], yr_data_lores[ii], zr_data_lores[ii], ramp, theta, phi, "\n", 
      file = con_lores, sep = ",")
}

# Write velocity receiver data
for (ii in 1:nrow(pos)) {
  cat(rtype2, xr_data_hires[ii], yr_data_hires[ii], zr_data_hires[ii], ramp, theta, phi, "\n", 
      file = con_hires, sep = ",")
  cat(rtype2, xr_data_lores[ii], yr_data_lores[ii], zr_data_lores[ii], ramp, theta, phi, "\n", 
      file = con_lores, sep = ",")
}

# Write acceleration receiver data
for (ii in 1:nrow(pos)) {
  cat(rtype3, xr_acc_hires[ii], yr_acc_hires[ii], zr_acc_hires[ii], ramp, theta, phi, "\n", 
      file = con_hires, sep = ",")
  cat(rtype3, xr_acc_lores[ii], yr_acc_lores[ii], zr_acc_lores[ii], ramp, theta, phi, "\n", 
      file = con_lores, sep = ",")
}

# Write pressure receiver data
for (ii in 1:nrow(pos)) {
  cat(rtype4, xr_data_hires[ii], yr_data_hires[ii], zr_data_hires[ii], ramp, theta, phi, "\n", 
      file = con_hires, sep = ",")
  cat(rtype4, xr_data_lores[ii], yr_data_lores[ii], zr_data_lores[ii], ramp, theta, phi, "\n", 
      file = con_lores, sep = ",")
}

# Plot output options
dscale = 1.0
tscale = 1.0
ascale = 1.0
iplot = 1

# Write output options
cat(dscale, tscale, ascale, iplot, "\n", file = con_hires, sep = ",")
cat(dscale, tscale, ascale, iplot, "\n", file = con_lores, sep = ",")

# Close file
close(con_hires)
close(con_lores)

