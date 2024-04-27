require("ggplot2")
require("animation")

#------------------------------------------------------
# Exact solution homogeneous linear MMS
#------------------------------------------------------
ExactSolHomoLin <- function(t, X) {
  K = 6.0e7
  G = 1.0e8
  alpha0 = 0.01
  mass = 1.0625000000000000e-07
  U = alpha0*X*t
  V = alpha0*X
  A = 0.0
  FF = 1 + alpha0*t
  sigma = (K + 4/3*G)*log(FF)
  B = 0
  return(data.frame(Time = t, Position = X, Disp = U, Vel = V, DefGrad = FF,
                    Stress = sigma, BodyForce = B))
}

ExactSolHomoQuad <- function(t, X) {
  K = 6.0e7
  G = 1.0e8
  alpha0 = 0.01
  mass = 1.0625000000000000e-07
  U = alpha0*X*t^2
  V = 2.0*alpha0*X*t
  A = 2.0*alpha0*X
  FF = 1 + alpha0*t^2
  sigma = (K + 4/3*G)*log(FF)
  B = 2.0*mass*alpha0*X
  return(data.frame(Time = t, Position = X, Disp = U, Vel = V, DefGrad = FF,
                    Stress = sigma, BodyForce = B))
}

ExactSolHarmonic <- function(t, X) {
  K = 6.0e7
  G = 1.0e8
  M = K + 4/3*G
  rho = 1700
  cp = sqrt(M/rho)
  alpha = 0.01
  omega = 1000
  mass = 1.0625000000000000e-07
  Urt = alpha*cos(omega*(t - X/cp))
  Uit = -alpha*sin(omega*(t - X/cp))
  U = Urt
  V = omega*Uit
  A = -omega^2*Urt
  FF = 1  - omega*Uit/cp
  sigma = (K + 4/3*G)*log(FF)
  B = mass*(omega^2*cp*Urt/(cp - omega*Uit) - omega^2*Urt)
  #sapply(t, function(time) {
  #  if (time == 0.0) {
  #    print(paste("t = ", time, " X = ", X, " U = ", U, " Uit = ", Uit , " Urt = ", Urt, " B = ", B, " cp = ", cp))
  #  }
  #})

  return(data.frame(Time = t, Position = X, Disp = U, Vel = V, DefGrad = FF,
                    Stress = sigma, BodyForce = B))
}

#------------------------------------------------------
# Extract the number of unique particles
#------------------------------------------------------
uniqueParticles <- function(data) {
  return(unique(data$id))
}

#------------------------------------------------------
# Extract a single curve
#------------------------------------------------------
extractCurve <- function(data, exact, initPos, t_index) {

  # Find the number of particles
  num.particles = length(uniqueParticles(data))
  start.index = (t_index-1)*num.particles+1
  end.index = t_index*num.particles
  print(paste("t_index", t_index, "startindex", start.index, "endindex", end.index))

  value = data$value[start.index:end.index]
  time = data$time[start.index:end.index]
  exact_value = exact$value[start.index:end.index]

  df1 = data.frame(Time = time, Position = initPos, Value = value, Label = "Simulation")
  df2 = data.frame(Time = time, Position = initPos, Value = exact_value, Label = "Exact")
  df = rbind(df1, df2)

  return(df)
}

#------------------------------------------------------
# Plot single curve
#------------------------------------------------------
plotCurve <- function(data, yLimits, dataType) {

  t_val = mean(data$Time) 
  tlabel = paste("Time = ", t_val*1.0e3, " ms")
  plt = ggplot(data = data) +
        geom_line(aes(x = Position, y = Value, color = Label), size = 1) +
        xlab("Initial Position") +
        ggtitle(tlabel) + 
        coord_cartesian(ylim = yLimits) + 
        theme_bw()

  if (dataType == "vel") {
    plt = plt + ylab("Particle velocity (m/s)") 
  } else if (dataType == "disp") {
    plt = plt + ylab("Particle displacement (m)") 
  } else if (dataType == "stress") {
    plt = plt + ylab("Particle stress (Pa)") 
  } else if (dataType == "defgrad") {
    plt = plt + ylab("Particle deformation gradient") 
  } else if (dataType == "fext") {
    plt = plt + ylab("External force (N)") 
  } else {
    print(paste("Unknown data type", dataType))
  }
  print(plt)
             
  #dev.copy(pdf, "gridVelBC.pdf")
  #dev.off()
}

#------------------------------------------------------
# Animate the plots
#------------------------------------------------------
animateCurve <- function(data, exact, yLimits, dataType) {

  # Find the number of particles
  num.particles = length(uniqueParticles(data))
  num.times = nrow(data)/num.particles

  # Find the initial positions of the particles
  initPos = data[1:num.particles,]$position

  lapply(seq(1, num.times, 1), function(t_index) {
    plotCurve(extractCurve(data, exact, initPos, t_index), yLimits, dataType)
  })
}

#------------------------------------------------------
# Function to read and plot the velocity/acc/stress data
#------------------------------------------------------
readAndPlotLargeDef <- function(inputDataFile, dispDataFile, ExactSol, yLimits, dataType) {

  outputGIFFile = paste0(unlist(strsplit(inputDataFile, "\\."))[1], ".gif")
  print(outputGIFFile)

  data = read.table(inputDataFile, header = FALSE, sep = "", 
                    stringsAsFactors=FALSE)

  dispdata = read.table(dispDataFile, header = FALSE, sep = "", 
                    stringsAsFactors=FALSE)

  if (dataType == "stress") {
    plotData = data.frame(time = data$V1, id = data$V4, 
                          value = data$V5, position = data$V14 - dispdata$V5) 
    exactVals = ExactSol(plotData$time, plotData$position)
    exactData = data.frame(time = exactVals$Time, position = exactVals$Position,
                           value = exactVals$Stress)
    print(paste("Stress Max/Min value = ", max(data$V5), ",", min(data$V5)))
  } else if (dataType == "defgrad") {
    plotData = data.frame(time = data$V1, id = data$V4, 
                          value = data$V5, position = data$V14 - dispdata$V5) 
    exactVals = ExactSol(plotData$time, plotData$position)
    exactData = data.frame(time = exactVals$Time, position = exactVals$Position,
                           value = exactVals$DefGrad)
    print(paste("Stress Max/Min value = ", max(data$V5), ",", min(data$V5)))
  } else if (dataType == "vel") {
    plotData = data.frame(time = data$V1, id = data$V4, 
                          value = data$V5, position = data$V8 - dispdata$V5) 
    exactVals = ExactSol(plotData$time, plotData$position)
    exactData = data.frame(time = exactVals$Time, position = exactVals$Position,
                           value = exactVals$Vel)
    print(paste("Vel Max/Min value = ", max(data$V5), ",", min(data$V5)))
  } else if (dataType == "fext") {
    plotData = data.frame(time = data$V1, id = data$V4, 
                          value = data$V5, position = data$V8 - dispdata$V5) 
    exactVals = ExactSol(plotData$time, plotData$position)
    exactData = data.frame(time = exactVals$Time, position = exactVals$Position,
                           value = exactVals$BodyForce)
    print(paste("Fext Max/Min value = ", max(data$V5), ",", min(data$V5)))
  } else if (dataType == "disp") {
    plotData = data.frame(time = data$V1, id = data$V4, 
                          value = data$V5, position = data$V8 - dispdata$V5) 
    exactVals = ExactSol(plotData$time, plotData$position)
    exactData = data.frame(time = exactVals$Time, position = exactVals$Position,
                           value = exactVals$Disp)
    print(paste("Disp Max/Min value = ", max(data$V5), ",", min(data$V5)))
  } else {
    print(paste("Unknown data type", dataType))
  }

  # Save as animation
  saveGIF(animateCurve(plotData, exactData, yLimits, dataType), interval=0.2, 
          movie.name=outputGIFFile)
}

#------------------------------------------------------
# Function to read and plot the velocity/acc/stress data
#------------------------------------------------------
readAndPlot <- function(inputDataFile, ExactSol, yLimits, dataType) {

  outputGIFFile = paste0(unlist(strsplit(inputDataFile, "\\."))[1], ".gif")
  print(outputGIFFile)

  data = read.table(inputDataFile, header = FALSE, sep = "", 
                    stringsAsFactors=FALSE)

  if (dataType == "stress") {
    plotData = data.frame(time = data$V1, id = data$V4, 
                          value = data$V5, position = data$V14) 
    exactVals = ExactSol(plotData$time, plotData$position)
    exactData = data.frame(time = exactVals$Time, position = exactVals$Position,
                           value = exactVals$Stress)
    print(paste("Stress Max/Min value = ", max(data$V5), ",", min(data$V5)))
  } else if (dataType == "defgrad") {
    plotData = data.frame(time = data$V1, id = data$V4, 
                          value = data$V5, position = data$V14) 
    exactVals = ExactSol(plotData$time, plotData$position)
    exactData = data.frame(time = exactVals$Time, position = exactVals$Position,
                           value = exactVals$DefGrad)
    print(paste("Stress Max/Min value = ", max(data$V5), ",", min(data$V5)))
  } else if (dataType == "vel") {
    plotData = data.frame(time = data$V1, id = data$V4, 
                          value = data$V5, position = data$V8) 
    exactVals = ExactSol(plotData$time, plotData$position)
    exactData = data.frame(time = exactVals$Time, position = exactVals$Position,
                           value = exactVals$Vel)
    print(paste("Vel Max/Min value = ", max(data$V5), ",", min(data$V5)))
  } else if (dataType == "fext") {
    plotData = data.frame(time = data$V1, id = data$V4, 
                          value = data$V5, position = data$V8) 
    exactVals = ExactSol(plotData$time, plotData$position)
    exactData = data.frame(time = exactVals$Time, position = exactVals$Position,
                           value = exactVals$BodyForce)
    print(paste("Fext Max/Min value = ", max(data$V5), ",", min(data$V5)))
  } else if (dataType == "disp") {
    plotData = data.frame(time = data$V1, id = data$V4, 
                          value = data$V5, position = data$V8) 
    exactVals = ExactSol(plotData$time, plotData$position)
    exactData = data.frame(time = exactVals$Time, position = exactVals$Position,
                           value = exactVals$Disp)
    print(paste("Disp Max/Min value = ", max(data$V5), ",", min(data$V5)))
  } else {
    print(paste("Unknown data type", dataType))
  }

  # Save as animation
  saveGIF(animateCurve(plotData, exactData, yLimits, dataType), interval=0.2, 
          movie.name=outputGIFFile)
}

#-------------------------------------------------------------------------
# Actually read and plot
#-------------------------------------------------------------------------
readAndPlotLargeDef("UniaxialStrainMMS_pDisplacement.dat", "UniaxialStrainMMS_pDisplacement.dat", 
             ExactSolHarmonic, c(-0.005, 0.011), "disp")
readAndPlotLargeDef("UniaxialStrainMMS_pF.dat", "UniaxialStrainMMS_pDisplacement.dat", 
             ExactSolHarmonic, c(0.9, 1.1), "defgrad")
readAndPlotLargeDef("UniaxialStrainMMS_pvel.dat", "UniaxialStrainMMS_pDisplacement.dat", 
             ExactSolHarmonic, c(-10, 10), "vel")
readAndPlotLargeDef("UniaxialStrainMMS_pfext.dat", "UniaxialStrainMMS_pDisplacement.dat", 
             ExactSolHarmonic, c(-5.0e-5, 5.0e-5), "fext")
readAndPlotLargeDef("UniaxialStrainMMS_pstress.dat", "UniaxialStrainMMS_pDisplacement.dat", 
             ExactSolHarmonic, c(-2.0e6, 6.0e6), "stress")

#readAndPlot("UniaxialStrainMMSHomoLin_pDisplacement.dat", 
#            ExactSolHomoLin, c(0.0, 1.1e-6), "disp")
#readAndPlot("UniaxialStrainMMSHomoLin_pF.dat", 
#            ExactSolHomoLin, c(0.9999, 1.0001), "defgrad")
#readAndPlot("UniaxialStrainMMSHomoLin_pvel.dat", 
#            ExactSolHomoLin, c(1.0e-6, 1.0e-3), "vel")
#readAndPlot("UniaxialStrainMMSHomoLin_pfext.dat", 
#            ExactSolHomoLin, c(-0.5e-5, 0.5e-5), "fext")
#readAndPlot("UniaxialStrainMMSHomoLin_pstress.dat", 
#            ExactSolHomoLin, c(0, 6000), "stress")

#readAndPlot("UniaxialStrainMMSHomoQuad_pDisplacement.dat", 
#            ExactSolHomoQuad, c(0.0, 1.1e-8), "disp")
#readAndPlot("UniaxialStrainMMSHomoQuad_pF.dat", 
#            ExactSolHomoQuad, c(0.9, 1.000010), "defgrad")
#readAndPlot("UniaxialStrainMMSHomoQuad_pvel.dat", 
#            ExactSolHomoQuad, c(0, 5.0e-6), "vel")
#readAndPlot("UniaxialStrainMMSHomoQuad_pfext.dat", 
#            ExactSolHomoQuad, c(-0.5e-9, 0.5e-9), "fext")
#readAndPlot("UniaxialStrainMMSHomoQuad_pstress.dat", 
#            ExactSolHomoQuad, c(0, 10), "stress")
