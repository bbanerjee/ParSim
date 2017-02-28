#--------------------------------------------------------------------------
# Name: convertPartDistToBinXML
#
# Purpose: Read the particle distribution csv file, and write out xml with 
#          ascii or base64 encoded data
#
# Usage:
#  source("convertPartDistToBinXML.R")
#  doConversion("ascii")
#  doConversion("base64")
#--------------------------------------------------------------------------

# Set the working directory and jre location
setwd(".")

# Clear existing variables
rm(list = ls())

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

# Install RCurl
if (!require(RCurl)) {
  install.packages("RCurl")
  library(RCurl)
}

# Function to increase indent size
xmlFormat <- function(xml, indent = 2) {
   s <- strsplit(saveXML(xml), "\n")[[1]]
   g <- gsubfn("^( +)", x ~ sprintf("%*s", indent * nchar(x), " "), s)
   paste(g, collapse = "\n")
}

doConversion <- function(partDistCSV, encoding = "ascii") {

  # Input particle distribution csv file
  #partDistCSV = "particle_distribution.csv"

  #------------------------------------------------------------------
  # Create XML
  # <?xml version='1.0' encoding='ISO-8859-1' ?>
  #<Ellip3D_input>
  #  <Particles number="100" type="ellipsoid">
  #    <id unit="none" encoding="base64" numComponents="1"> 1 2 3 ... 100 </id>
  #    <radii unit="m" encoding="base64" numComponents="3" > a1 b1 c1 a2 b2 c2 ... </radii>
  #    <position unit="m" encoding="base64" numComponents="3" > x1 y1 z1 x2 y2 z2 ... </position>
  #    <axle_a_angle unit="radian" encoding="base64" numComponents="3" > theta_x1 theta_y1 theta_z1 ... </axle_a_angle>
  #    <axle_b_angle unit="radian" encoding="base64" numComponents="3" > theta_x1 theta_y1 theta_z1 ... </axle_b_angle>
  #    <axle_c_angle unit="radian" encoding="base64" numComponents="3"> theta_x1 theta_y1 theta_z1 ... </axle_c_angle>
  #    <velocity unit="m/s" encoding="base64" numComponents="3"> x1 y1 z1 x2 y2 z2 ... </velocity>
  #    <omega unit="radian/s" encoding="base64" numComponents="3"> x1 y1 z1 x2 y2 z2 ... </omega>
  #    <force unit="N" encoding="base64" numComponents="3"> x1 y1 z1 x2 y2 z2 ... </force>
  #    <moment unit="Nm" encoding="base64" numComponents="3"> x1 y1 z1 x2 y2 z2 ... </moment>
  #  </Particles>
  #  <Sieves number="5">
  #    <percent_passing> 0.1 0.2 0.3 ... <percent_passing>
  #    <size unit="mm"> 0.1 0.2 0.3 ... <size>
  #    <sieve_ratio>
  #      <ratio_ba>0.8</ratio_ba>
  #      <ratio_ca>0.6</ratio_ca>
  #    </sieve_ratio>
  #  </Sieves>
  #</Ellip3D_input>
  #------------------------------------------------------------------

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

  xml <- xmlOutputDOM(tag = "Ellip3D_input")

  partTypes = unique(df_part$type)
  for (partType in partTypes) {
    partTypeName = "ellipsoid"
    if (partType == 1) {
      partTypeName = "sphere"
    } else if (partType == 2) {
      partTypeName = "cube"
    } else if (partType == 3) {
      partTypeName = "polyellipsoid"
    } else {
      partTypeName = "ellipsoid"
    }
    df = df_part[which(df_part$type == partType),]
    xml$addTag("Particles", attrs = c(number = nrow(df), type = partTypeName), close = FALSE)

    # Add Id
    partID = paste(df$id, collapse=" ")
    if (encoding == "base64") {
      partID = base64(partID)
    }
    xml$addTag("id", partID, attrs = c(unit = "none", encoding = encoding, numComponents = "1"))

    # Add radius
    radius = cbind(df$radius_a, df$radius_b, df$radius_c)
    partVec = paste(t(radius), collapse=" ")
    if (encoding == "base64") {
      partVec = base64(partVec)
    }
    xml$addTag("radii", partVec, attrs = c(unit = "m", encoding = encoding, numComponents = "3"))

    # Add axle_a
    axle_a = cbind(df$axle_a_x, df$axle_a_y, df$axle_a_z)
    partVec = paste(t(axle_a), collapse=" ")
    if (encoding == "base64") {
      partVec = base64(partVec)
    }
    xml$addTag("axle_a", partVec, attrs = c(unit = "rad", encoding = encoding, numComponents = "3"))

    # Add axle_b
    axle_b = cbind(df$axle_b_x, df$axle_b_y, df$axle_b_z)
    partVec = paste(t(axle_b), collapse=" ")
    if (encoding == "base64") {
      partVec = base64(partVec)
    }
    xml$addTag("axle_b", partVec, attrs = c(unit = "rad", encoding = encoding, numComponents = "3"))

    # Add axle_c
    axle_c = cbind(df$axle_c_x, df$axle_c_y, df$axle_c_z)
    partVec = paste(t(axle_c), collapse=" ")
    if (encoding == "base64") {
      partVec = base64(partVec)
    }
    xml$addTag("axle_c", partVec, attrs = c(unit = "rad", encoding = encoding, numComponents = "3"))

    # Add position
    position = cbind(df$position_x, df$position_y, df$position_z)
    partVec = paste(t(position), collapse=" ")
    if (encoding == "base64") {
      partVec = base64(partVec)
    }
    xml$addTag("position", partVec, attrs = c(unit = "m", encoding = encoding, numComponents = "3"))

    # Add velocity
    velocity = cbind(df$velocity_x, df$velocity_y, df$velocity_z)
    partVec = paste(t(velocity), collapse=" ")
    if (encoding == "base64") {
      partVec = base64(partVec)
    }
    xml$addTag("velocity", partVec, attrs = c(unit = "m/s", encoding = encoding, numComponents = "3"))

    # Add omega
    omega = cbind(df$omga_x, df$omga_y, df$omga_z)
    partVec = paste(t(omega), collapse=" ")
    if (encoding == "base64") {
      partVec = base64(partVec)
    }
    xml$addTag("omega", partVec, attrs = c(unit = "rad/s", encoding = encoding, numComponents = "3"))

    # Add force
    force = cbind(df$force_x, df$force_y, df$force_z)
    partVec = paste(t(force), collapse=" ")
    if (encoding == "base64") {
      partVec = base64(partVec)
    }
    xml$addTag("force", partVec, attrs = c(unit = "N", encoding = encoding, numComponents = "3"))

    # Add moment
    moment = cbind(df$moment_x, df$moment_y, df$moment_z)
    partVec = paste(t(moment), collapse=" ")
    if (encoding == "base64") {
      partVec = base64(partVec)
    }
    xml$addTag("moment", partVec, attrs = c(unit = "Nm", encoding = encoding, numComponents = "3"))

    xml$closeTag() # Close the Particles tage
  }

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
  percent = paste(df_sieve$V1, collapse = " ")
  size = paste(df_sieve$V2, collapse = " ")
  xml$addTag("percent_passing", percent)
  xml$addTag("size", size, attrs = c(unit = "mm"))
  xml$addTag("sieve_ratio", close = FALSE)
  xml$addTag("ratio_ba", df_sieve_ratio[1])
  xml$addTag("ratio_ca", df_sieve_ratio[2])
  xml$closeTag() # sieve_ratio
  xml$closeTag() # Sieves

  #------------------------------------------------------------------
  # Save the XML file
  #------------------------------------------------------------------
  partDistXML = ""
  if (encoding == "ascii") {
    partDistXML = paste0(sub(".csv", "", partDistCSV), ".ascii.xml")
  } else {
    partDistXML = paste0(sub(".csv", "", partDistCSV), ".base64.xml")
  }
  cat(xmlFormat(xml), sep = "\n", file = partDistXML)
}

doConversion("particle_distribution.short.csv", "ascii")
doConversion("particle_distribution.short.csv", "base64")
doConversion("particle_distribution.csv", "ascii")
doConversion("particle_distribution.csv", "base64")


