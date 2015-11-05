require("ggplot2")
require("animation")

#------------------------------------------------------
# Extract the number of cell indices 
#------------------------------------------------------
extractCellIndices <- function(data) {
  return(unique(data$V1))
}

#------------------------------------------------------
# Extract a single curve
#------------------------------------------------------
extractCurve <- function(velBCData, t_index) {

  cell.indices = extractCellIndices(velBCData)
  num.cell.indices = length(cell.indices)
  start.index = (t_index-1)*num.cell.indices+1
  end.index = t_index*num.cell.indices
  print(paste("t_index", t_index, "startindex", start.index, "endindex", end.index))

  v_velBC = velBCData$V4[start.index:end.index]
  cell = velBCData$V1[start.index:end.index]

  return(data.frame(Cell = cell.indices, VelBC = v_velBC))
}

#------------------------------------------------------
# Plot single curve
#------------------------------------------------------
plotCurve <- function(data, yLimits, dataType) {

  
  plt = ggplot(data = data) +
        geom_line(aes(x = Cell, y = VelBC), color = "red", size = 1) +
        xlab("Node index") +
        coord_cartesian(ylim = yLimits) + 
        theme_bw()

  if (dataType == "pvel") {
    plt = plt + ylab("Particle velocity (m/s)") 
  } else if (dataType == "pstress") {
    plt = plt + ylab("Particle stress (Pa)") 
  } else if (dataType == "gacc") {
    plt = plt + ylab("Grid acceleration (m/s^2)") 
  } else {
    plt = plt + ylab("Grid velocity (m/s)") 
  }
  print(plt)
             
  #dev.copy(pdf, "gridVelBC.pdf")
  #dev.off()
}

#------------------------------------------------------
# Animate the plots
#------------------------------------------------------
animateCurve <- function(velBCData, yLimits, dataType) {
  num.cell.indices = length(extractCellIndices(velBCData))
  num.times = nrow(velBCData)/num.cell.indices
  lapply(seq(1, num.times, 1), function(t_index) {
    plotCurve(extractCurve(velBCData, t_index), yLimits, dataType)
  })
}

#------------------------------------------------------
# Function to read and plot the velocity/acc/stress data
#------------------------------------------------------
readAndPlot <- function(inputDataFile, yLimits, dataType) {

  outputGIFFile = paste0(unlist(strsplit(inputDataFile, "\\."))[1], ".gif")
  print(outputGIFFile)

  if (dataType == "pstress") {
    v_VelBC = read.table(inputDataFile, skip = 2, header = FALSE, sep = "", 
                         stringsAsFactors=FALSE)

    print(paste("Stress Max/Min value = ", max(v_VelBC$V4), ",", min(v_VelBC$V4)))
  } else {
    v_VelBC = read.table(inputDataFile, skip = 2, header = FALSE, sep = "", 
                         stringsAsFactors=FALSE)
    v_VelBC$V4 = gsub("\\[","",v_VelBC$V4)
    v_VelBC$V6 = gsub("\\]","",v_VelBC$V6)
    v_VelBC$V4 = as.numeric(v_VelBC$V4)
    v_VelBC$V6 = as.numeric(v_VelBC$V6)

    print(paste("Vel/Acc Max/Min value = ", max(v_VelBC$V4), ",", min(v_VelBC$V4)))
  }

  # Save as animation
  saveGIF(animateCurve(v_VelBC, yLimits, dataType), interval=0.2, 
          movie.name=outputGIFFile)
}

#-------------------------------------------------------------------------
# Actually read and plot
#-------------------------------------------------------------------------
#readAndPlot("undamped_midres_tractionBC_gacc.dat",
#               c(-200, 200), "gacc")
#readAndPlot("undamped_midres_tractionBC_gvel.dat",
#               c(-0.003, 0.003), "gvel")
#readAndPlot("undamped_midres_tractionBC_pvel.dat",
#               c(-0.003, 0.003), "pvel")
#readAndPlot("undamped_midres_tractionBC_pstress.dat",
#               c(-500, 500), "pstress")
#readAndPlot("undamped_midres_velBC_pstress.dat",
#            c(-1000000, 0), "pstress")
#readAndPlot("undamped_momform_midres_velBC_gvel.dat",
#            c(-0.5, 1.5), "gvel")
#readAndPlot("undamped_midres_impact_velBC_pvel.dat",
#            c(-0.5, 1.5), "pvel")
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#readAndPlot("damped_midres_impact_velBC_gvel.dat",
#            c(-0.5, 1.5), "gvel")
#readAndPlot("damped_midres_impact_velBC_pvel.dat",
#            c(-0.5, 1.5), "pvel")
#readAndPlot("damped_midres_tractionBC_pstress.dat",
#            c(-500, 500), "pstress")
#readAndPlot("damped_midres_velBC_gvel.dat",
#            c(-0.5, 1.5), "gvel")
#readAndPlot("damped_midres_velBC_pvel.dat",
#            c(-0.5, 1.5), "pvel")
#readAndPlot("damped_midres_velBC_pstress.dat",
#            c(-1000000, 0), "pstress")
#readAndPlot("damped_uintah_midres_velBC_gvel.dat",
#            c(-0.5, 1.5), "gvel")

#readAndPlot("damped_midres_SquareVelBC_pvel.dat", 
#               c(-0.1, 1.1), "pvel")
#readAndPlot("damped_midres_SquareVelBC.dat", 
#               c(-0.1, 1.1), "gvel")
#readAndPlot("damped_midres_SquarePressFromVelBC.dat", 
#               c(-0.1, 1.1), "gvel")
#readAndPlot("damped_midres_HatVelBC.dat", 
#               c(-0.1, 1), "gvel")
#readAndPlot("damped_midres_PressFromVelBC.dat", 
#               c(-0.1, 1), "gvel")
#readAndPlot("undamped_arenisca_lores_PressFromVelBC.dat", 
#               c(-0.25, 1), "gvel")
#readAndPlot("damped_arenisca_lores_PressFromVelBC.dat", 
#               c(-0.25, 1), "gvel")
#readAndPlot("undamped_arenisca_midres_PressFromVelBC.dat", 
#               c(-0.1, 1), "gvel")
#readAndPlot("damped_arenisca_midres_PressFromVelBC.dat", 
#               c(-0.1, 1), "gvel")
#readAndPlot("undamped_lores_explosionBC.dat", 
#               c(-50, 250), "gvel")
#readAndPlot("damped_lores_explosionBC.dat", 
#               c(-50, 250), "gvel")
#readAndPlot("undamped_midres_explosionBC.dat", 
#               c(-50, 250), "gvel")
#readAndPlot("damped_midres_explosionBC.dat",
#               c(-50, 250), "gvel")
#readAndPlot("damped_hires_explosionBC.dat",
#               c(-50, 250), "gvel")
#readAndPlot("damped_lores_cfl_explosionBC.dat",
#               c(-50, 250), "gvel")
#readAndPlot("damped_lores_cpdi_explosionBC.dat",
#               c(-50, 250), "gvel")


