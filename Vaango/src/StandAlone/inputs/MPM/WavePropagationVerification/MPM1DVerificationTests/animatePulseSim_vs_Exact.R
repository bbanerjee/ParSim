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
ExactSol <- function(t) {
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
extractCurve <- function(velBCData, pressBCData, t_index) {

  v = 1
  dt = 1.0e-5
  t = t_index*dt
  dx = 0.001
  x_exact = ExactSol(t)
  x_index = round(x_exact/dx)
  
  cell.indices = extractCellIndices(pressBCData)
  num.cell.indices = length(cell.indices)
  start.index = (t_index-1)*num.cell.indices+1
  end.index = t_index*num.cell.indices
  print(paste("t_index", t_index, "startindex", start.index, "endindex", end.index))

  v_pressBC = pressBCData$V4[start.index:end.index]
  cell = pressBCData$V1[start.index:end.index]

  v_exact = v_pressBC
  v_exact[1:x_index] = v
  v_exact[(x_index+1):num.cell.indices] = 0

  v_velBC = velBCData$V4

  return(data.frame(Cell = cell.indices, Exact = v_exact, VelBC = v_velBC, PressBC = v_pressBC))
}

#------------------------------------------------------
# Plot single curve
#------------------------------------------------------
plotCurve <- function(data) {

  plt = ggplot(data = data) +
        geom_line(aes(x = Cell, y = Exact), size = 1) +
        geom_line(aes(x = Cell, y = VelBC), color = "red", size = 1) +
        geom_line(aes(x = Cell, y = PressBC), color = "blue", size = 1) +
        xlab("Node index") +
        ylab("Grid velocity (m/s)") +
        coord_cartesian(ylim = c(-0.1, 1.3)) + 
        theme_bw()
  print(plt)
             
  #dev.copy(pdf, "gridVelBC.pdf")
  #dev.off()
}

#------------------------------------------------------
# Animate the plots
#------------------------------------------------------
animateCurve <- function(velBCData, pressBCData) {
  num.cell.indices = length(extractCellIndices(pressBCData))
  num.times = nrow(pressBCData)/num.cell.indices
  lapply(seq(1, num.times, 1), function(t_index) {
    plotCurve(extractCurve(velBCData, pressBCData, t_index))
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
saveGIF(animateCurve(v_VelBC, v_PressBC), interval=0.1, movie.name="damped_midres_square_compare.gif")

