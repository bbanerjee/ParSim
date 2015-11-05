#--------------------------------------------------------------------------
# Name: writeReactionLoadCurve
#
# Purpose: Read the recation curve produced by vaaango, 
#          write out xml and plot load as a function of time
#
# Usage:
#  source("writeReactionLoadCurve.R")
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

# Load curve csv file
loadCurveCSV = "RigidReactionForceHat.dat"
loadCurveXML = "RigidReactionForceHat.xml"
loadCurvePDF = "RigidReactionForceHat.pdf"
loadCurvePNG = "RigidReactionForceHat.png"

# Area over which load is applied
len = 2*0.0005
area = len*len

# df
df = read.table(loadCurveCSV, header = FALSE, sep = "", stringsAsFactors=FALSE)
df$V2 = gsub("\\[","",df$V2)
df$V4 = gsub("\\]","",df$V4)
df$V2 = as.numeric(df$V2)
df$V4 = as.numeric(df$V4)

# Compute pressure
df$V2 = df$V2/area

# Create XML
# <?xml version='1.0' encoding='ISO-8859-1' ?>
#<Uintah_Include>
#          <load_curve>
#            <id>1</id>
#            <time_point>
#              <time> 0 </time>
#              <load> 0 </load>
#            </time_point>
#xml <- xmlOutputDOM(tag = "Uintah_Include")
#xml$addTag("load_curve", close = FALSE)
#xml$addTag("id", 1)
#for (ii in 1:nrow(df)) {
#  xml$addTag("time_point", close = FALSE)
#  xml$addTag("time", df[ii,1])
#  xml$addTag("load", df[ii,2])
#  xml$closeTag()
#}
#xml$closeTag()

xml <- xmlTree(tag = "Uintah_Include")
xml$addNode("load_curve", close = FALSE)
xml$addNode("id", 1)
for (ii in 1:nrow(df)) {
  xml$addNode("time_point", close = FALSE)
  xml$addNode("time", df[ii,1])
  xml$addNode("load", df[ii,2])
  xml$closeTag()
}
xml$closeTag()

# Save the XML file
saveXML(xml$value(), file = loadCurveXML, prefix = '<?xml version="1.0" encoding="ISO-8859-1" ?>\n')

# Read the load curve file file
ups_data <- xmlParse(loadCurveXML)

# Convert the data using xpath
time_path = "//time_point/time"
load_path = "//time_point/load"
df <- data.frame(
  time = as.numeric(xpathSApply(ups_data, time_path, xmlValue)),
  load = as.numeric(xpathSApply(ups_data, load_path, xmlValue)))

# Plot the load curve
yticks = pretty(log10(df$load*1.0e-6), 10)
yticks = yticks[which(round(yticks) == yticks)]

plt = ggplot(data = df, aes(x = time*1.0e3, y = load*1.0e-6)) +
      geom_line(color = "#0072B2", size = 1) +
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
dev.copy(png, loadCurvePNG)
dev.off()



