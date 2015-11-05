# Set the working directory
setwd(".")

# Load required packages
require("ggplot2")

# Function to read stress data
readStressData <- function(stressFileName) {
  aa = read.table(stressFileName, skip = 2, header = FALSE, sep = "", stringsAsFactors=FALSE)
  #plot(aa$V4)
  #lines(aa$V4, col=2)
  return(aa)
}

# Function to read velocity data
readVelocityData <- function(velocityFileName) {
  aa = read.table(velocityFileName, skip = 2, header = FALSE, sep = "", stringsAsFactors=FALSE)
  aa$V4 = gsub("\\[","",aa$V4)
  aa$V6 = gsub("\\]","",aa$V6)
  aa$V4 = as.numeric(aa$V4)
  aa$V6 = as.numeric(aa$V6)
  return(aa)
}

# Exact solution for an input step function
Exact <- function(t) {
  rho = 1700
  K = 6.0e7
  G = 9.0e7
  M = K + 4.0/3.0*G
  c_p = sqrt(M/rho)
  x0 = 0
  x = x0 + c_p*t
  return(x)
}

# Read and plot
fileName1 = "OneD_impact_velocityBC_gimp_damped.dat"
fileName2 = "OneD_pressureBC_from_velBC_gimp_damped.uda.dat"
aa = readVelocityData(fileName1)  # Generated using a velocity BC
bb = readVelocityData(fileName2)  # Generated using a pressure BC derived from the
                                  # RigidReactionForce file generated in the 
                                  # velocity BC run

v = 1
t = 0.00100001
x_exact = Exact(t)
dx = 0.001
index = round(x_exact/dx)
aa$v_exact = aa$V4
aa$v_exact[1:index] = v
aa$v_exact[(index+1):nrow(aa)] = 0
aa$v_press = bb$V4

df1 = data.frame(x = aa$V1, y = aa$V4, label = "Velocity BC")
df2 = data.frame(x = aa$V1, y = aa$v_press, label = "Pressure BC from Reaction")
df3 = data.frame(x = aa$V1, y = aa$v_exact, label = "Exact solution")
df = rbind(df1, df2, df3)
plt = ggplot(data = df) +
      geom_line(aes(x = x, y = y, color = label)) + 
      xlab("Node index") +
      ylab("Grid velocity (m/s)") +
      theme_bw()
print(plt)
             
#plot(aa$V4)
#lines(aa$V4, col=2)
dev.copy(pdf, "comparisonSimExactGridVelBC.pdf")
dev.off()
