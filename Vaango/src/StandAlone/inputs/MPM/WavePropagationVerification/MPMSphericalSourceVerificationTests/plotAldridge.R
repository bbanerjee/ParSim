# Set wd
setwd(".")

# Clean up
rm(list = ls())
for (i in dev.list()) {dev.off(i)}

# Load ggplot2
require("ggplot2")

# load colorspace
require("colorspace")

simAccFile_1 = "../Spherical_source_verification_gimp.uda.000_acc.data"
simDataFile_1 = "../Spherical_source_verification_gimp.uda.000.partdata"
aldridgeOutFile_1 = "output_data_lores.out"
runCase_1 = "Full GIMP Lores"

simAccFile_2 = "../Spherical_source_verification_gimp_quarter.uda.000_acc.data"
simDataFile_2 = "../Spherical_source_verification_gimp_quarter.uda.000.partdata"
aldridgeOutFile_2 = "output_data_lores.out"
runCase_2 = "Eighth GIMP Lores"

simAccFile_3 = "../Spherical_source_verification_gimp_quarter_hires.uda.000_acc.data"
simDataFile_3 = "../Spherical_source_verification_gimp_quarter_hires.uda.000.partdata"
aldridgeOutFile_3 = "output_data_hires.out"
runCase_3 = "Eighth GIMP Hires"

simAccFile_4 = "../Spherical_source_verification_gimp_quarter_hires_cpdi.uda.000_acc.data"
simDataFile_4 = "../Spherical_source_verification_gimp_quarter_hires_cpdi.uda.000.partdata"
aldridgeOutFile_4 = "output_data_hires.out"
runCase_4 = "Eighth CPDI CBDI Hires"

simAccFile_5 = "../Spherical_source_verification_gimp_quarter_hires_cpdi_nocbdi.uda.000_acc.data"
simDataFile_5 = "../Spherical_source_verification_gimp_quarter_hires_cpdi_nocbdi.uda.000.partdata"
aldridgeOutFile_5 = "output_data_hires.out"
runCase_5 = "Eighth CPDI no CBDI Hires"

simAccFile_6 = "../Spherical_source_verification_gimp_quarter_hires_cpdi_nocbdi_nodamp.uda.000_acc.data"
simDataFile_6 = "../Spherical_source_verification_gimp_quarter_hires_cpdi_nocbdi_nodamp.uda.000.partdata"
aldridgeOutFile_6 = "output_data_hires.out"
runCase_6 = "Eighth CPDI no CBDI no Damping Hires"

#---------------------------------------------------------------------------
# Read exact data
#---------------------------------------------------------------------------
readAldridgeData <- function(aldridgeOutFile, runCase) {

  # Read the output file generated by explode_for.exe
  tt = read.csv(aldridgeOutFile, header = FALSE, sep = "", comment.char=">")
  tt$RunCase = runCase

  # The data contains the accelerations and the pressures at each of the sensor
  # locations
  names(tt) = c("Time", "Position", "Value", "RunCase")

  # If the positions vary then we cannot use position as a factor
  # Just extract the data by row index
  nrows = nrow(tt)
  exact_subsets = split(tt, rep(1:24, each=nrows/24))
  
  # Extract disp, vel, accelerations and pressures
  sensor1_disp  = exact_subsets[[1]]
  sensor2_disp  = exact_subsets[[2]]
  sensor3_disp  = exact_subsets[[3]]
  sensor4_disp  = exact_subsets[[4]]
  sensor5_disp  = exact_subsets[[5]]
  sensor6_disp  = exact_subsets[[6]]
  sensor1_vel  = exact_subsets[[7]]
  sensor2_vel  = exact_subsets[[8]]
  sensor3_vel  = exact_subsets[[9]]
  sensor4_vel  = exact_subsets[[10]]
  sensor5_vel  = exact_subsets[[11]]
  sensor6_vel  = exact_subsets[[12]]
  sensor1_acc  = exact_subsets[[13]]
  sensor2_acc  = exact_subsets[[14]]
  sensor3_acc  = exact_subsets[[15]]
  sensor4_acc  = exact_subsets[[16]]
  sensor5_acc  = exact_subsets[[17]]
  sensor6_acc  = exact_subsets[[18]]
  sensor1_press  = exact_subsets[[19]]
  sensor2_press  = exact_subsets[[20]]
  sensor3_press  = exact_subsets[[21]]
  sensor4_press  = exact_subsets[[22]]
  sensor5_press  = exact_subsets[[23]]
  sensor6_press  = exact_subsets[[24]]

  # Create a list to return
  data_list = list(sensor1_disp = sensor1_disp,
                   sensor2_disp = sensor2_disp,
                   sensor3_disp = sensor3_disp,
                   sensor4_disp = sensor4_disp,
                   sensor5_disp = sensor5_disp,
                   sensor6_disp = sensor6_disp,
                   sensor1_vel = sensor1_vel,
                   sensor2_vel = sensor2_vel,
                   sensor3_vel = sensor3_vel,
                   sensor4_vel = sensor4_vel,
                   sensor5_vel = sensor5_vel,
                   sensor6_vel = sensor6_vel,
                   sensor1_acc = sensor1_acc,
                   sensor2_acc = sensor2_acc,
                   sensor3_acc = sensor3_acc,
                   sensor4_acc = sensor4_acc,
                   sensor5_acc = sensor5_acc,
                   sensor6_acc = sensor6_acc,
                   sensor1_press = sensor1_press,
                   sensor2_press = sensor2_press,
                   sensor3_press = sensor3_press,
                   sensor4_press = sensor4_press,
                   sensor5_press = sensor5_press,
                   sensor6_press = sensor6_press)
  
  # Colors
  #mycolors = rainbow_hcl(6)
  
  # Displacement plot
  #dev.new()
  #gg1 = ggplot() + 
  #      geom_line(data = sensor1_disp, aes(x=Time, y=Value), color = mycolors[1]) + 
  #      geom_line(data = sensor2_disp, aes(x=Time, y=Value), color = mycolors[2]) + 
  #      geom_line(data = sensor3_disp, aes(x=Time, y=Value), color = mycolors[3]) + 
  #      geom_line(data = sensor4_disp, aes(x=Time, y=Value), color = mycolors[4]) + 
  #      geom_line(data = sensor5_disp, aes(x=Time, y=Value), color = mycolors[5]) + 
  #      geom_line(data = sensor6_disp, aes(x=Time, y=Value), color = mycolors[6]) + 
  #      xlim(0,0.002)
  #print(gg1)
  
  # Velocity plot
  #dev.new()
  #gg2 = ggplot() + 
  #      geom_line(data = sensor1_vel, aes(x=Time, y=Value), color = mycolors[1]) + 
  #      geom_line(data = sensor2_vel, aes(x=Time, y=Value), color = mycolors[2]) + 
  #      geom_line(data = sensor3_vel, aes(x=Time, y=Value), color = mycolors[3]) + 
  #      geom_line(data = sensor4_vel, aes(x=Time, y=Value), color = mycolors[4]) + 
  #      geom_line(data = sensor5_vel, aes(x=Time, y=Value), color = mycolors[5]) + 
  #      geom_line(data = sensor6_vel, aes(x=Time, y=Value), color = mycolors[6]) + 
  #      xlim(0,0.002)
  #print(gg2)
  
  # Acceleration plot
  #dev.new()
  #gg3 = ggplot() + 
  #      geom_line(data = sensor1_acc, aes(x=Time, y=Value), color = mycolors[1]) + 
  #      geom_line(data = sensor2_acc, aes(x=Time, y=Value), color = mycolors[2]) + 
  #      geom_line(data = sensor3_acc, aes(x=Time, y=Value), color = mycolors[3]) + 
  #      geom_line(data = sensor4_acc, aes(x=Time, y=Value), color = mycolors[4]) + 
  #      geom_line(data = sensor5_acc, aes(x=Time, y=Value), color = mycolors[5]) + 
  #      geom_line(data = sensor6_acc, aes(x=Time, y=Value), color = mycolors[6]) + 
  #      xlim(0,0.002)
  #print(gg3)
  
  # Pressure plot
  #dev.new()
  #gg4 = ggplot() + 
  #      geom_line(data = sensor1_press, aes(x=Time, y=Value), color = mycolors[1]) + 
  #      geom_line(data = sensor2_press, aes(x=Time, y=Value), color = mycolors[2]) + 
  #      geom_line(data = sensor3_press, aes(x=Time, y=Value), color = mycolors[3]) + 
  #      geom_line(data = sensor4_press, aes(x=Time, y=Value), color = mycolors[4]) + 
  #      geom_line(data = sensor5_press, aes(x=Time, y=Value), color = mycolors[5]) + 
  #      geom_line(data = sensor6_press, aes(x=Time, y=Value), color = mycolors[6]) + 
  #      xlim(0,0.002)
  #print(gg4)
  
  return(data_list)
}

#---------------------------------------------------------------------------
# Rbind data from two run cases
#---------------------------------------------------------------------------
rbindAldridgeData <- function(case1, case2) {
  data_list = list(sensor1_disp = rbind(case1$sensor1_disp, case2$sensor1_disp),
                   sensor2_disp = rbind(case1$sensor2_disp, case2$sensor2_disp),
                   sensor3_disp = rbind(case1$sensor3_disp, case2$sensor3_disp),
                   sensor4_disp = rbind(case1$sensor4_disp, case2$sensor4_disp),
                   sensor5_disp = rbind(case1$sensor5_disp, case2$sensor5_disp),
                   sensor6_disp = rbind(case1$sensor6_disp, case2$sensor6_disp),
                   sensor1_vel = rbind(case1$sensor1_vel, case2$sensor1_vel),
                   sensor2_vel = rbind(case1$sensor2_vel, case2$sensor2_vel),
                   sensor3_vel = rbind(case1$sensor3_vel, case2$sensor3_vel),
                   sensor4_vel = rbind(case1$sensor4_vel, case2$sensor4_vel),
                   sensor5_vel = rbind(case1$sensor5_vel, case2$sensor5_vel),
                   sensor6_vel = rbind(case1$sensor6_vel, case2$sensor6_vel),
                   sensor1_acc = rbind(case1$sensor1_acc, case2$sensor1_acc),
                   sensor2_acc = rbind(case1$sensor2_acc, case2$sensor2_acc),
                   sensor3_acc = rbind(case1$sensor3_acc, case2$sensor3_acc),
                   sensor4_acc = rbind(case1$sensor4_acc, case2$sensor4_acc),
                   sensor5_acc = rbind(case1$sensor5_acc, case2$sensor5_acc),
                   sensor6_acc = rbind(case1$sensor6_acc, case2$sensor6_acc),
                   sensor1_press = rbind(case1$sensor1_press, case2$sensor1_press),
                   sensor2_press = rbind(case1$sensor2_press, case2$sensor2_press),
                   sensor3_press = rbind(case1$sensor3_press, case2$sensor3_press),
                   sensor4_press = rbind(case1$sensor4_press, case2$sensor4_press),
                   sensor5_press = rbind(case1$sensor5_press, case2$sensor5_press),
                   sensor6_press = rbind(case1$sensor6_press, case2$sensor6_press))
  return(data_list)
}

#---------------------------------------------------------------------------
# Read the MPM simulation data
#---------------------------------------------------------------------------
readMPMData <- function(simAccFile, simDataFile, runCase) {

  # Read the simulated accelerations
  acc_data = read.csv(simAccFile, header = FALSE, sep = "")

  # Create a subset of columns
  acc_data = acc_data[,5:9]

  # Add run case
  acc_data$RunCase = runCase
  names(acc_data) = c("Time", "xacc", "yacc", "zacc", "Position", "RunCase")

  # Convert position into a factor
  acc_data$Position = as.factor(acc_data$Position)
  
  # Create subsets based on position
  subsets = split(acc_data, acc_data$Position, drop = TRUE)
  sim_sensor1_acc = subsets[[1]]
  sim_sensor2_acc = subsets[[2]]
  sim_sensor3_acc = subsets[[3]]
  sim_sensor4_acc = subsets[[4]]
  sim_sensor5_acc = subsets[[5]]
  sim_sensor6_acc = subsets[[6]]
  
  # Plot the accelerations
  #dev.new()
  #gg3 = gg3 +
  #      geom_line(data = sim_sensor1_acc, aes(x=Time, y=xacc), color = mycolors[1]) + 
  #      geom_line(data = sim_sensor2_acc, aes(x=Time, y=xacc), color = mycolors[2]) + 
  #      geom_line(data = sim_sensor3_acc, aes(x=Time, y=xacc), color = mycolors[3]) + 
  #      geom_line(data = sim_sensor4_acc, aes(x=Time, y=xacc), color = mycolors[4]) + 
  #      geom_line(data = sim_sensor5_acc, aes(x=Time, y=xacc), color = mycolors[5]) +
  #      geom_line(data = sim_sensor6_acc, aes(x=Time, y=xacc), color = mycolors[6])  
  #print(gg3)

  # Read the simulated positions, velocities, and stresses
  part_data = read.csv(simDataFile, header = FALSE, sep = "")
  
  # Add run case
  part_data$RunCase = runCase
  names(part_data) = c("Time", "x", "y", "z", "vx", "vy", "vz",  "sigx", "sigy", "sigz", "RunCase")
  
  # Compute mean stress
  part_data$press = (part_data$sigx + part_data$sigy + part_data$sigz)/3.0
  
  # Extract disp, vel, and stresses
  nrows = nrow(part_data)
  nrow_sensor = nrows/6
  start_1 = 1
  end_1 = nrow_sensor
  start_2 = end_1+1
  end_2 = 2*nrow_sensor
  start_3 = end_2+1
  end_3 = 3*nrow_sensor
  start_4 = end_3+1
  end_4 = 4*nrow_sensor
  start_5 = end_4+1
  end_5 = 5*nrow_sensor
  start_6 = end_5+1
  end_6 = nrows
  
  sim_sensor1_data  = part_data[start_1:end_1,]
  sim_sensor2_data  = part_data[start_2:end_2,]
  sim_sensor3_data  = part_data[start_3:end_3,]
  sim_sensor4_data  = part_data[start_4:end_4,]
  sim_sensor5_data  = part_data[start_5:end_5,]
  sim_sensor6_data  = part_data[start_6:end_6,]
  
  # Compute displacement
  sim_sensor1_data$ux = sim_sensor1_data$x - sim_sensor1_data$x[1]
  sim_sensor1_data$uy = sim_sensor1_data$y - sim_sensor1_data$y[1]
  sim_sensor1_data$uz = sim_sensor1_data$z - sim_sensor1_data$z[1]
  sim_sensor2_data$ux = sim_sensor2_data$x - sim_sensor2_data$x[1]
  sim_sensor2_data$uy = sim_sensor2_data$y - sim_sensor2_data$y[1]
  sim_sensor2_data$uz = sim_sensor2_data$z - sim_sensor2_data$z[1]
  sim_sensor3_data$ux = sim_sensor3_data$x - sim_sensor3_data$x[1]
  sim_sensor3_data$uy = sim_sensor3_data$y - sim_sensor3_data$y[1]
  sim_sensor3_data$uz = sim_sensor3_data$z - sim_sensor3_data$z[1]
  sim_sensor4_data$ux = sim_sensor4_data$x - sim_sensor4_data$x[1]
  sim_sensor4_data$uy = sim_sensor4_data$y - sim_sensor4_data$y[1]
  sim_sensor4_data$uz = sim_sensor4_data$z - sim_sensor4_data$z[1]
  sim_sensor5_data$ux = sim_sensor5_data$x - sim_sensor5_data$x[1]
  sim_sensor5_data$uy = sim_sensor5_data$y - sim_sensor5_data$y[1]
  sim_sensor5_data$uz = sim_sensor5_data$z - sim_sensor5_data$z[1]
  sim_sensor6_data$ux = sim_sensor6_data$x - sim_sensor6_data$x[1]
  sim_sensor6_data$uy = sim_sensor6_data$y - sim_sensor6_data$y[1]
  sim_sensor6_data$uz = sim_sensor6_data$z - sim_sensor6_data$z[1]
  
  # Create a list to return
  data_list = list(sensor1_acc = sim_sensor1_acc,
                   sensor2_acc = sim_sensor2_acc,
                   sensor3_acc = sim_sensor3_acc,
                   sensor4_acc = sim_sensor4_acc,
                   sensor5_acc = sim_sensor5_acc,
                   sensor6_acc = sim_sensor6_acc,
                   sensor1_data = sim_sensor1_data,
                   sensor2_data = sim_sensor2_data,
                   sensor3_data = sim_sensor3_data,
                   sensor4_data = sim_sensor4_data,
                   sensor5_data = sim_sensor5_data,
                   sensor6_data = sim_sensor6_data)
  # Plot
  #gg1 = gg1 +
  #      geom_line(data = sim_sensor1_data, aes(x=Time, y=ux), color = mycolors[1]) + 
  #      geom_line(data = sim_sensor2_data, aes(x=Time, y=ux), color = mycolors[2]) + 
  #      geom_line(data = sim_sensor3_data, aes(x=Time, y=ux), color = mycolors[3]) + 
  #      geom_line(data = sim_sensor4_data, aes(x=Time, y=ux), color = mycolors[4]) + 
  #      geom_line(data = sim_sensor5_data, aes(x=Time, y=ux), color = mycolors[5]) + 
  #      geom_line(data = sim_sensor6_data, aes(x=Time, y=ux), color = mycolors[6]) 
  #print(gg1)
  
  #gg2 = gg2 +
  #      geom_line(data = sim_sensor1_data, aes(x=Time, y=vx), color = mycolors[1]) + 
  #      geom_line(data = sim_sensor2_data, aes(x=Time, y=vx), color = mycolors[2]) + 
  #      geom_line(data = sim_sensor3_data, aes(x=Time, y=vx), color = mycolors[3]) + 
  #      geom_line(data = sim_sensor4_data, aes(x=Time, y=vx), color = mycolors[4]) + 
  #      geom_line(data = sim_sensor5_data, aes(x=Time, y=vx), color = mycolors[5]) + 
  #      geom_line(data = sim_sensor6_data, aes(x=Time, y=vx), color = mycolors[6]) 
  #print(gg2)
  
  #gg4 = gg4 +
  #      geom_line(data = sim_sensor1_data, aes(x=Time, y=press), color = mycolors[1]) + 
  #      geom_line(data = sim_sensor2_data, aes(x=Time, y=press), color = mycolors[2]) + 
  #      geom_line(data = sim_sensor3_data, aes(x=Time, y=press), color = mycolors[3]) + 
  #      geom_line(data = sim_sensor4_data, aes(x=Time, y=press), color = mycolors[4]) + 
  #      geom_line(data = sim_sensor5_data, aes(x=Time, y=press), color = mycolors[5]) + 
  #      geom_line(data = sim_sensor6_data, aes(x=Time, y=press), color = mycolors[6]) 
  #print(gg4)

  return(data_list)
}

#---------------------------------------------------------------------------
# Rbind data from two run cases
#---------------------------------------------------------------------------
rbindMPMData <- function(case1, case2) {
  data_list = list(sensor1_acc = rbind(case1$sensor1_acc, case2$sensor1_acc),
                   sensor2_acc = rbind(case1$sensor2_acc, case2$sensor2_acc),
                   sensor3_acc = rbind(case1$sensor3_acc, case2$sensor3_acc),
                   sensor4_acc = rbind(case1$sensor4_acc, case2$sensor4_acc),
                   sensor5_acc = rbind(case1$sensor5_acc, case2$sensor5_acc),
                   sensor6_acc = rbind(case1$sensor6_acc, case2$sensor6_acc),
                   sensor1_data = rbind(case1$sensor1_data, case2$sensor1_data),
                   sensor2_data = rbind(case1$sensor2_data, case2$sensor2_data),
                   sensor3_data = rbind(case1$sensor3_data, case2$sensor3_data),
                   sensor4_data = rbind(case1$sensor4_data, case2$sensor4_data),
                   sensor5_data = rbind(case1$sensor5_data, case2$sensor5_data),
                   sensor6_data = rbind(case1$sensor6_data, case2$sensor6_data))
  return(data_list)
}

#-----------------------------------------------
# Function for plotting the displacement
#-----------------------------------------------
plotDisp <- function(exact_disp, sim_disp, sensor_id, outFile) {
  
  exact_df = data.frame(Time = exact_disp$Time*1.0e3, ux = exact_disp$Value*1.0e3,
                  Label = paste("Aldridge", exact_disp$RunCase))
  sim_df = data.frame(Time = sim_disp$Time*1.0e3, ux = sim_disp$ux*1.0e3,
                  Label = paste("MPM", sim_disp$RunCase))
  df = rbind(exact_df, sim_df)
  
  dev.new()
  plt = ggplot(data = df) + 
        geom_line(aes(x=Time, y=ux, color=Label), size=1) + 
        xlim(0,2) +
        ggtitle(paste("Sensor ", sensor_id, " radius = ", sim_disp$x[1])) +
        xlab("Time (ms)") +
        ylab("Displacement (mm)") +
        scale_color_discrete(name = "") +
        theme_bw() +
        theme(legend.justification=c(1,0), legend.position=c(1,0),
            plot.title = element_text(size = 10),
            axis.title.x = element_text(size=16),
            axis.title.y = element_text(size=16),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14),
            legend.text = element_text(size=12))
  print(plt)
  dev.copy(pdf, outFile)
  dev.off()

  return(plt)
}

#-----------------------------------------------
# Function for plotting the velocity
#-----------------------------------------------
plotVel <- function(exact_vel, sim_vel, sensor_id, outFile) {
  
  exact_df = data.frame(Time = exact_vel$Time*1.0e3, vx = exact_vel$Value*1.0e3,
                  Label = paste("Aldridge", exact_vel$RunCase))
  sim_df = data.frame(Time = sim_vel$Time*1.0e3, vx = sim_vel$vx*1.0e3,
                  Label = paste("MPM", sim_vel$RunCase))
  df = rbind(exact_df, sim_df)
  
  dev.new()
  plt = ggplot(data = df) + 
        geom_line(aes(x=Time, y=vx, color=Label), size=1) + 
        xlim(0,2) +
        ggtitle(paste("Sensor ", sensor_id, " radius = ", sim_vel$x[1])) +
        xlab("Time (ms)") +
        ylab("Velocity (m/s)") +
        scale_color_discrete(name = "") +
        theme_bw() +
        theme(legend.justification=c(1,0), legend.position=c(1,0),
            plot.title = element_text(size = 10),
            axis.title.x = element_text(size=16),
            axis.title.y = element_text(size=16),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14),
            legend.text = element_text(size=12))
  print(plt)
  dev.copy(pdf, outFile)
  dev.off()

  return(plt)
}

#-----------------------------------------------
# Function for plotting the acceleration
#-----------------------------------------------
plotAcc <- function(exact_acc, sim_acc, sensor_id, outFile) {

  print(head(exact_acc))
  print(head(sim_acc))
  exact_df = data.frame(Time = exact_acc$Time*1.0e3, acc = exact_acc$Value/9.81,
                  Label = paste("Aldridge", exact_acc$RunCase))
  sim_df = data.frame(Time = sim_acc$Time*1.0e3, acc = sim_acc$xacc/9.81,
                  Label = paste("MPM", sim_acc$RunCase))
  df = rbind(exact_df, sim_df)
  print(summary(df))
  
  dev.new()
  plt = ggplot(data = df) + 
        geom_line(aes(x=Time, y=acc, color=Label), size=1) + 
        xlim(0,2) +
        ggtitle(paste("Sensor ", sensor_id, " radius = ", sim_acc$Position[1])) +
        xlab("Time (ms)") +
        ylab("Acceleration (g)") +
        scale_color_discrete(name = "") +
        theme_bw() +
        theme(legend.justification=c(1,0), legend.position=c(1,0),
            plot.title = element_text(size = 10),
            axis.title.x = element_text(size=16),
            axis.title.y = element_text(size=16),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14),
            legend.text = element_text(size=12))
  print(plt)
  dev.copy(pdf, outFile)
  dev.off()

  return(plt)
}

#-----------------------------------------------
# Function for plotting the pressure
#-----------------------------------------------
plotPress <- function(exact_press, sim_press, sensor_id, outFile) {

  exact_df = data.frame(Time = exact_press$Time*1.0e3, press = exact_press$Value*1.0e-6,
                  Label = paste("Aldridge", exact_press$RunCase))
  sim_df = data.frame(Time = sim_press$Time*1.0e3, press = sim_press$press*1.0e-6,
                  Label = paste("MPM", sim_press$RunCase))
  df = rbind(exact_df, sim_df)
  
  dev.new()
  plt = ggplot(data = df) + 
        geom_line(aes(x=Time, y=press, color=Label), size=1) + 
        xlim(0,2) +
        ggtitle(paste("Sensor ", sensor_id, " radius = ", sim_press$x[1])) +
        xlab("Time (ms)") +
        ylab("Pressure (MPa)") +
        scale_color_discrete(name = "") +
        theme_bw() +
        theme(legend.justification=c(1,0), legend.position=c(1,0),
            plot.title = element_text(size = 10),
            axis.title.x = element_text(size=16),
            axis.title.y = element_text(size=16),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14),
            legend.text = element_text(size=12))
  print(plt)
  dev.copy(pdf, outFile)
  dev.off()

  return(plt)
}

#---------------------------------------------------------------------------
# Plot the data
#---------------------------------------------------------------------------
plotComparisons <- function(exact, sim, simAccFile, simDataFile) {

  outFile = paste0(simAccFile, "_1.pdf")
  plotAcc(exact$sensor1_acc, mpm$sensor1_acc, 1, outFile)
  outFile = paste0(simAccFile, "_2.pdf")
  plotAcc(exact$sensor2_acc, mpm$sensor2_acc, 2, outFile)
  outFile = paste0(simAccFile, "_3.pdf")
  plotAcc(exact$sensor3_acc, mpm$sensor3_acc, 3, outFile)
  outFile = paste0(simAccFile, "_4.pdf")
  plotAcc(exact$sensor4_acc, mpm$sensor4_acc, 4, outFile)
  outFile = paste0(simAccFile, "_5.pdf")
  plotAcc(exact$sensor5_acc, mpm$sensor5_acc, 5, outFile)
  outFile = paste0(simAccFile, "_6.pdf")
  plotAcc(exact$sensor6_acc, mpm$sensor6_acc, 6, outFile)
  
  #---------------------------------------------------------------------------
  outFile = paste0(simDataFile, "disp_1.pdf")
  plotDisp(exact$sensor1_disp, mpm$sensor1_data, 1, outFile)
  outFile = paste0(simDataFile, "disp_2.pdf")
  plotDisp(exact$sensor2_disp, mpm$sensor2_data, 2, outFile)
  outFile = paste0(simDataFile, "disp_3.pdf")
  plotDisp(exact$sensor3_disp, mpm$sensor3_data, 3, outFile)
  outFile = paste0(simDataFile, "disp_4.pdf")
  plotDisp(exact$sensor4_disp, mpm$sensor4_data, 4, outFile)
  outFile = paste0(simDataFile, "disp_5.pdf")
  plotDisp(exact$sensor5_disp, mpm$sensor5_data, 5, outFile)
  outFile = paste0(simDataFile, "disp_6.pdf")
  plotDisp(exact$sensor6_disp, mpm$sensor6_data, 6, outFile)
  
  outFile = paste0(simDataFile, "vel_1.pdf")
  plotVel(exact$sensor1_vel, mpm$sensor1_data, 1, outFile)
  outFile = paste0(simDataFile, "vel_2.pdf")
  plotVel(exact$sensor2_vel, mpm$sensor2_data, 2, outFile)
  outFile = paste0(simDataFile, "vel_3.pdf")
  plotVel(exact$sensor3_vel, mpm$sensor3_data, 3, outFile)
  outFile = paste0(simDataFile, "vel_4.pdf")
  plotVel(exact$sensor4_vel, mpm$sensor4_data, 4, outFile)
  outFile = paste0(simDataFile, "vel_5.pdf")
  plotVel(exact$sensor5_vel, mpm$sensor5_data, 5, outFile)
  outFile = paste0(simDataFile, "vel_6.pdf")
  plotVel(exact$sensor6_vel, mpm$sensor6_data, 6, outFile)
  
  outFile = paste0(simDataFile, "press_1.pdf")
  plotPress(exact$sensor1_press, mpm$sensor1_data, 1, outFile)
  outFile = paste0(simDataFile, "press_2.pdf")
  plotPress(exact$sensor2_press, mpm$sensor2_data, 2, outFile)
  outFile = paste0(simDataFile, "press_3.pdf")
  plotPress(exact$sensor3_press, mpm$sensor3_data, 3, outFile)
  outFile = paste0(simDataFile, "press_4.pdf")
  plotPress(exact$sensor4_press, mpm$sensor4_data, 4, outFile)
  outFile = paste0(simDataFile, "press_5.pdf")
  plotPress(exact$sensor5_press, mpm$sensor5_data, 5, outFile)
  outFile = paste0(simDataFile, "press_6.pdf")
  plotPress(exact$sensor6_press, mpm$sensor6_data, 6, outFile)
}

#---------------------------------------------------------------------------
# Process  data
#---------------------------------------------------------------------------
exact_1 = readAldridgeData(aldridgeOutFile_1, runCase_1)
mpm_1 = readMPMData(simAccFile_1, simDataFile_1, runCase_1)
exact_2 = readAldridgeData(aldridgeOutFile_2, runCase_2)
mpm_2 = readMPMData(simAccFile_2, simDataFile_2, runCase_2)

exact = rbindAldridgeData(exact_1, exact_2)
mpm = rbindMPMData(mpm_1, mpm_2)

plotComparisons(exact, mpm, simAccFile_2, simDataFile_2)

#exact_1 = readAldridgeData(aldridgeOutFile_2, runCase_2)
#mpm_1 = readMPMData(simAccFile_2, simDataFile_2, runCase_2)
#exact_2 = readAldridgeData(aldridgeOutFile_3, runCase_3)
#mpm_2 = readMPMData(simAccFile_3, simDataFile_3, runCase_3)

#exact = rbindAldridgeData(exact_1, exact_2)
#mpm = rbindMPMData(mpm_1, mpm_2)

#plotComparisons(exact, mpm, simAccFile_3, simDataFile_3)

#exact_1 = readAldridgeData(aldridgeOutFile_3, runCase_3)
#mpm_1 = readMPMData(simAccFile_3, simDataFile_3, runCase_3)
#exact_2 = readAldridgeData(aldridgeOutFile_4, runCase_4)
#mpm_2 = readMPMData(simAccFile_4, simDataFile_4, runCase_4)

#exact = rbindAldridgeData(exact_1, exact_2)
#mpm = rbindMPMData(mpm_1, mpm_2)

#plotComparisons(exact, mpm, simAccFile_4, simDataFile_4)

#exact_1 = readAldridgeData(aldridgeOutFile_4, runCase_4)
#mpm_1 = readMPMData(simAccFile_4, simDataFile_4, runCase_4)
#exact_2 = readAldridgeData(aldridgeOutFile_5, runCase_5)
#mpm_2 = readMPMData(simAccFile_5, simDataFile_5, runCase_5)

#exact = rbindAldridgeData(exact_1, exact_2)
#mpm = rbindMPMData(mpm_1, mpm_2)

#plotComparisons(exact, mpm, simAccFile_5, simDataFile_5)

#exact_1 = readAldridgeData(aldridgeOutFile_5, runCase_5)
#mpm_1 = readMPMData(simAccFile_5, simDataFile_5, runCase_5)
#exact_2 = readAldridgeData(aldridgeOutFile_6, runCase_6)
#mpm_2 = readMPMData(simAccFile_6, simDataFile_6, runCase_6)

#exact = rbindAldridgeData(exact_1, exact_2)
#mpm = rbindMPMData(mpm_1, mpm_2)

#plotComparisons(exact, mpm, simAccFile_6, simDataFile_6)

