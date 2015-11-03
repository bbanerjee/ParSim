#--------------------------------------------------------------------------
# Name: plotLoadCurve
#
# Purpose: Read the UPS load curve xml and plot load as a function of time
#
# Usage:
#  source("plotLoadCurve.R")
#--------------------------------------------------------------------------

# Set the working directory and jre location
setwd(".")
#Sys.setenv(JAVA_HOME="/usr/lib/jvm/java-7-openjdk-amd64/jre")

# Install rjava
if (!require(rJava)) {
  install.packages("rJava")
  library(rJava)
}

# Install XML
if (!require(XML)) {
  install.packages("XML")
  library(XML)
}

# Install ggplot2 
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library("ggplot2")
}

# Install scales 
if (!require("scales")) {
  install.packages("scales")
  library("scales")
}

# Install data.table
if (!require(data.table)) {
  install.packages("data.table")
  library(data.table)
}

# Install pracma  (for interp1)
if (!require(pracma)) {
  install.packages("pracma")
  library("pracma")
}

# Load curve file
#loadCurveFile = "loadCurve.xml"
#loadCurveFile = "loadCurveDelay.xml"
#loadCurveFile = "LoadCurveNoBucketInitStress.xml"
#loadCurveFile = "LoadCurveDelayOffset.xml"
#loadCurveFile = "CentrifugeLoadCurveHigh.xml"
loadCurveFile = "CentrifugeLoadCurveLow.xml"
#loadCurveFile = "CentrifugeLoadCurveMid.xml"

# Read the load curve file file and replace data
ups_data <- xmlParse(loadCurveFile)

# Convert the data using xpath
time_path = "//time_point/time"
load_path = "//time_point/load"
df <- data.frame(
  time = as.numeric(xpathSApply(ups_data, time_path, xmlValue)),
  load = as.numeric(xpathSApply(ups_data, load_path, xmlValue)))

# Plot the load curve
yticks = pretty(log10(df$load*1.0e-6), 10)
yticks = yticks[which(round(yticks) == yticks)]

plt = ggplot() + 
      geom_line(data = df, aes(x = time*1.0e3, y = load*1.0e-6), color = "#0072B2", size = 1) +
      geom_point(color = "#D55E00") +
      xlab("Time (millisec)") + 
      ylab("Load (MPa)") + 
      theme_bw() + 
      xlim(0, 10) +
      scale_y_log10(limits = c(1.0-3, 1.0e3),
                    breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) + 
      theme(legend.justification=c(1,0), legend.position=c(1,0),
            plot.title = element_text(size = 10),
            axis.title.x = element_text(size=16),
            axis.title.y = element_text(size=16),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14),
            legend.text = element_text(size=12))
print(plt)
#dev.copy(png, "load_curve_hi.png")
#dev.copy(png, "load_curve_delay.png")
#dev.copy(png, "load_curve_nobucket_initstress.png")
#dev.copy(png, "load_curve_delay_offset.png")
#dev.copy(png, "load_curve_high.png")
dev.copy(png, "load_curve_low.png")
#dev.copy(png, "load_curve_mid.png")
dev.off()

# Write the load curve in Aldridge format
# after resampling 
minTime = min(df$time)
maxTime = max(df$time)
numSamples = 40000
timeSample = seq(minTime, maxTime, length.out = numSamples)
loadSample = interp1(df$time, df$load, timeSample, method="linear")
resampledDF = data.frame(time = timeSample, load = loadSample)
outFile = paste(loadCurveFile, "resampled", sep=".")
write.table(resampledDF, file = outFile, sep = " ",
            row.names = FALSE, col.names = FALSE)

dev.new()
plt1 = ggplot() + 
      geom_line(data = df, 
                aes(x = time*1.0e3, y = load*1.0e-6), color = "#A172B2", size = 1) +
      geom_line(data = resampledDF, 
                aes(x = time*1.0e3, y = load*1.0e-6), color = "#0072B2", size = 1) +
      geom_point(color = "#D55E00") +
      xlab("Time (millisec)") + 
      ylab("Load (MPa)") + 
      theme_bw() + 
      xlim(0, 10) +
      scale_y_log10(limits = c(1.0-3, 1.0e3),
                    breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) + 
      theme(legend.justification=c(1,0), legend.position=c(1,0),
            plot.title = element_text(size = 10),
            axis.title.x = element_text(size=16),
            axis.title.y = element_text(size=16),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14),
            legend.text = element_text(size=12))
print(plt1)




