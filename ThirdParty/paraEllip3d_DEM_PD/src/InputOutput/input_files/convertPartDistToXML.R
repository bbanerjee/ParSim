#--------------------------------------------------------------------------
# Name: convertPartDistToXML
#
# Purpose: Read the particle distribution csv file, and write out xml
#
# Usage:
#  source("convertPartDistToXML.R")
#--------------------------------------------------------------------------

# Set the working directory and jre location
setwd(".")

# Install rjava
#if (!require(rJava)) {
#  install.packages("rJava")
#  library(rJava)
#}

# Install XML
if (!require(XML)) {
  install.packages("XML")
  library(XML)
}

# Install gsubfn
if (!require(gsubfn)) {
  install.packages("gsubfn")
  library(gsubfn)
}

# Function to increase indent size
xmlFormat <- function(xml, indent = 2) {
   s <- strsplit(saveXML(xml), "\n")[[1]]
   g <- gsubfn("^( +)", x ~ sprintf("%*s", indent * nchar(x), " "), s)
   paste(g, collapse = "\n")
}

# Input particle distribution csv file
partDistCSV = "particle_distribution.csv"

#------------------------------------------------------------------
# First handle the particle data
#------------------------------------------------------------------
# Read the first line (contains the number of particles)
df_numPart = read.csv(partDistCSV, header = FALSE, nrows = 1, sep="", 
                   blank.lines.skip = TRUE)
numPart = df_numPart$V1

# Read the particle data into a data frame
# Data should contain the following header:
# > names(df)
# [1] "id"         "type"       "radius_a"   "radius_b"   "radius_c"  
# [6] "position_x" "position_y" "position_z" "axle_a_x"   "axle_a_y"  
#[11] "axle_a_z"   "axle_b_x"   "axle_b_y"   "axle_b_z"   "axle_c_x"  
#[16] "axle_c_y"   "axle_c_z"   "velocity_x" "velocity_y" "velocity_z"
#[21] "omga_x"     "omga_y"     "omga_z"     "force_x"    "force_y"   
#[26] "force_z"    "moment_x"   "moment_y"   "moment_z"  
df_part = read.csv(partDistCSV, header = TRUE, skip = 1, sep="", 
                   blank.lines.skip = TRUE, nrows = numPart)

# Create XML
# <?xml version='1.0' encoding='ISO-8859-1' ?>
#<Ellip3D_input>
#  <Particles number="100">
#    <particle type="ellipsoid" id="1">
#      <radii> [a b c] </radii>
#      <position> [x y z] </position>
#      <axle_a_angle> [theta_x theta_y theta_z] </axle_a_angle>
#      <axle_b_angle> [theta_x theta_y theta_z] </axle_b_angle>
#      <axle_c_angle> [theta_x theta_y theta_z] </axle_c_angle>
#      <velocity> [x y z] </velocity>
#      <omega> [x y z] </omega>
#      <force> [x y z] </force>
#      <moment> [x y z] </moment>
#    </particle>
#  </Particles>
#</Ellip3D_input>
xml <- xmlOutputDOM(tag = "Ellip3D_input")
xml$addTag("Particles", attrs = c(number = numPart), close = FALSE)
for (ii in 1:nrow(df_part)) {
  partID = df$id[ii];
  partTypeNum = df$type[ii];
  partType = ""
  if (partTypeNum == 0) {
    partType = "ellipsoid"
  }
  xml$addTag("particle", attrs = c(type = partType, id = partID),  close = FALSE)

  # Add radius
  partVec = paste("[", df$radius_a[ii], ",", df$radius_b[ii], ",", df$radius_c[ii], "]");
  xml$addTag("radii", partRadius, attrs = c(unit = "m"))

  # Add axle_a
  partVec = paste("[", df$axle_a_x[ii], ",", df$axle_a_y[ii], ",", df$axle_a_z[ii], "]");
  xml$addTag("axle_a_angle", partVec, attrs = c(unit = "rad"))

  # Add axle_b
  partVec = paste("[", df$axle_b_x[ii], ",", df$axle_b_y[ii], ",", df$axle_b_z[ii], "]");
  xml$addTag("axle_b_angle", partVec, attrs = c(unit = "rad"))

  # Add axle_c
  partVec = paste("[", df$axle_c_x[ii], ",", df$axle_c_y[ii], ",", df$axle_c_z[ii], "]");
  xml$addTag("axle_c_angle", partVec, attrs = c(unit = "rad"))

  # Add position
  partVec = paste("[", df$position_x[ii], ",", df$position_y[ii], ",", df$position_z[ii], "]");
  xml$addTag("position", partVec, attrs = c(unit = "m"))

  # Add velocity
  partVec = paste("[", df$velocity_x[ii], ",", df$velocity_y[ii], ",", df$velocity_z[ii], "]");
  xml$addTag("velocity", partVec, attrs = c(unit = "m/s"))

  # Add omega
  partVec = paste("[", df$omga_x[ii], ",", df$omga_y[ii], ",", df$omga_z[ii], "]");
  xml$addTag("omega", partVec, attrs = c(unit = "rad/s"))

  # Add force
  partVec = paste("[", df$force_x[ii], ",", df$force_y[ii], ",", df$force_z[ii], "]");
  xml$addTag("force", partVec, attrs = c(unit = "N"))

  # Add moment
  partVec = paste("[", df$moment_x[ii], ",", df$moment_y[ii], ",", df$moment_z[ii], "]");
  xml$addTag("moment", partVec, attrs = c(unit = "Nm"))

  xml$closeTag()
}
xml$closeTag()

#------------------------------------------------------------------
# Next handle the sieve data
#------------------------------------------------------------------
# Read the number of sieves
df_numSieve = read.csv(partDistCSV, header = FALSE, skip = numPart+2, sep="", 
                       blank.lines.skip = TRUE, nrows = 1)
numSieve = df_numSieve$V1
# Read the sieves
df_sieve = read.csv(partDistCSV, header = FALSE, skip = (numPart+4), sep="", 
                    blank.lines.skip = TRUE, nrows = numSieve)
# Read the sieve ratios
df_sieve_ratio = read.csv(partDistCSV, header = FALSE, skip = numPart+numSieve+4, sep="", 
                          blank.lines.skip = TRUE, nrows = 1)

# Create sieve xml
xml$addTag("Sieves", attrs = c(number = numSieve), close = FALSE)
for (ii in 1:nrow(df_sieve)) {
  xml$addTag("sieve", attrs = c(id = ii), close = FALSE)
  xml$addTag("percent_passing", df_sieve[[1]][ii])
  xml$addTag("size", df_sieve[[2]][ii], attrs = c(unit = "mm"))
  xml$closeTag();
}
{
  xml$addTag("sieve_ratio", close = FALSE)
  xml$addTag("ratio_ba", df_sieve_ratio[1])
  xml$addTag("ratio_ca", df_sieve_ratio[2])
  xml$closeTag
}
xml$closeTag()

#------------------------------------------------------------------
# Save the XML file
#------------------------------------------------------------------
partDistXML = paste0(sub(".csv", "", partDistCSV), ".xml")
cat(xmlFormat(xml), sep = "\n", file = partDistXML, append = TRUE)

