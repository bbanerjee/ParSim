require("ggplot2")
require("animation")

#------------------------------------------------------
# Extract the number of cell indices 
#------------------------------------------------------
extractCellIndices <- function(data) {
  return(unique(data$V1))
}

#------------------------------------------------------
# Compute the exact solution for a given time
#------------------------------------------------------
ExactSolHat <- function(t) {
  rho = 1700
  K = 6.0e7
  G = 9.0e7
  M = K + 4.0/3.0*G
  c_p = sqrt(M/rho)
  x0 = 0
  x = x0 + c_p*t
  return(x)
}

#------------------------------------------------------
# Extract a single curve
#------------------------------------------------------
extractCurve <- function(velBCData, pressBCData, exactFunction, t_index) {

  cell.indices = extractCellIndices(pressBCData)
  num.cell.indices = length(cell.indices)
  start.index = (t_index-1)*num.cell.indices+1
  end.index = t_index*num.cell.indices
  print(paste("t_index", t_index, "startindex", start.index, "endindex", end.index))

  cell = pressBCData$V1[start.index:end.index]
  v_pressBC = pressBCData$V4[start.index:end.index]
  v_velBC = velBCData$V4[start.index:end.index]

  v_exact = exactFunction(t_index, num.cell.indices)

  df_velBC = data.frame(Cell = cell.indices, Velocity = v_velBC, Label = 'Velocity BC')
  df_pressBC = data.frame(Cell = cell.indices, Velocity = v_pressBC, Label = 'Pressure BC')
  df_exact = data.frame(Cell = cell.indices, Velocity = v_exact, Label = 'Exact')
  df_all = rbind(df_velBC, df_pressBC, df_exact);

  return(df_all)
}

#------------------------------------------------------
# Compute the exact solution for a given time
#------------------------------------------------------
ExactSolVelHat <- function(t_index, num.pts) {

  # Material properties
  rho = 1700
  K = 6.0e7
  G = 9.0e7
  M = K + 4.0/3.0*G
  c_p = sqrt(M/rho)

  # Discretization
  dt = 1.0e-5
  dx = 0.001

  # Current time
  t = t_index*dt

  # Velocity pulse
  velocityFunction <- function(v0, x) {

    vel = x
    vel[which(x < 0)] = v0
    vel[which(x >= 0)] = 0
    return(vel)
  }

  # Input x values
  xx = seq(0, (num.pts-1)*dx, dx)

  # Current position
  xx = xx - c_p*t

  # Compute velocity
  v0 = 1
  v_exact = velocityFunction(v0, xx)

  return(v_exact)
}

ExactSolVelTri <- function(t_index, num.pts) {

  # Velocity pulse
  velocityPulse <- function(tt) {

    t_pts = c(-1.0, 0, 0.0001, 0.0002, 1.0)
    v_pts = c(0.0, 0.0, 1.0, 0.0, 0.0) 
  
    ind0 = tail(which(t_pts <= tt), n=1)
    ind1 = head(which(t_pts > tt), n=1)
    t0 = t_pts[ind0]
    t1 = t_pts[ind1]
    v0 = v_pts[ind0]
    v1 = v_pts[ind1]

    s = (tt - t0)/(t1 - t0)
    v = (1-s)*v0 + s*v1
    return(v)
  }
    
  # Material properties
  rho = 1700
  K = 6.0e7
  G = 9.0e7
  M = K + 4.0/3.0*G
  c_p = sqrt(M/rho)

  # Discretization
  dt = 1.0e-5
  dx = 0.001

  # Input x values
  xx = seq(0.0, (num.pts-1)*dx, dx)

  # Current time
  t = t_index*dt

  # Compute velocity
  t_shift = t - xx/c_p
  v_exact = sapply(t_shift, function(tt) {velocityPulse(tt)})

  return(v_exact)
}

#------------------------------------------------------
# Plot single curve
#------------------------------------------------------
plotCurve <- function(data) {

  plt = ggplot(data = data) +
        geom_path(aes(x = Cell, y = Velocity, color = Label), size = 1) +
        #geom_line(aes(x = Cell, y = VelBC), color = "red", size = 1) +
        #geom_line(aes(x = Cell, y = PressBC), color = "blue", size = 1) +
        xlab("Node index") +
        ylab("Grid velocity (m/s)") +
        coord_cartesian(ylim = c(-0.1, 1.3)) + 
        theme_bw() +
        theme(legend.position="top")
  print(plt)
             
  #dev.copy(pdf, "gridVelBC.pdf")
  #dev.off()
}

#------------------------------------------------------
# Animate the plots
#------------------------------------------------------
animateCurve <- function(velBCData, pressBCData, exactFunction) {
  num.cell.indices = length(extractCellIndices(pressBCData))
  num.times = nrow(pressBCData)/num.cell.indices
  lapply(seq(1, num.times, 1), function(t_index) {
    plotCurve(extractCurve(velBCData, pressBCData, exactFunction, t_index))
  })
}

#------------------------------------------------------
# Read the data
#------------------------------------------------------
v_VelBC = read.table("damped_midres_SquareVelBC.dat", skip = 2, header = FALSE, sep = "", stringsAsFactors=FALSE)
v_VelBC$V4 = gsub("\\[","",v_VelBC$V4)
v_VelBC$V6 = gsub("\\]","",v_VelBC$V6)
v_VelBC$V4 = as.numeric(v_VelBC$V4)
v_VelBC$V6 = as.numeric(v_VelBC$V6)

v_PressBC = read.table("damped_midres_SquarePressFromVelBC.dat", skip = 2, header = FALSE, sep = "", stringsAsFactors=FALSE)
v_PressBC$V4 = gsub("\\[","",v_PressBC$V4)
v_PressBC$V6 = gsub("\\]","",v_PressBC$V6)
v_PressBC$V4 = as.numeric(v_PressBC$V4)
v_PressBC$V6 = as.numeric(v_PressBC$V6)

# Save as animation
saveGIF(animateCurve(v_VelBC, v_PressBC, ExactSolVelHat), interval=0.1, movie.name="damped_midres_square_compare.gif")


v_VelBC = read.table("damped_midres_HatVelBC.dat", skip = 2, header = FALSE, sep = "", stringsAsFactors=FALSE)
v_VelBC$V4 = gsub("\\[","",v_VelBC$V4)
v_VelBC$V6 = gsub("\\]","",v_VelBC$V6)
v_VelBC$V4 = as.numeric(v_VelBC$V4)
v_VelBC$V6 = as.numeric(v_VelBC$V6)

v_PressBC = read.table("damped_midres_PressFromVelBC.dat", skip = 2, header = FALSE, sep = "", stringsAsFactors=FALSE)
v_PressBC$V4 = gsub("\\[","",v_PressBC$V4)
v_PressBC$V6 = gsub("\\]","",v_PressBC$V6)
v_PressBC$V4 = as.numeric(v_PressBC$V4)
v_PressBC$V6 = as.numeric(v_PressBC$V6)

saveGIF(animateCurve(v_VelBC, v_PressBC, ExactSolVelTri), interval=0.1, movie.name="damped_midres_hat_compare.gif")

