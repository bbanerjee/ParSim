#--------------------------------------------------------------------------
# Name: convertPartDistToBinXML
#
# Purpose: Read the boundary csv file, and write out xml with 
#          ascii or base64 encoded data
#
# Usage:
#  source("convertPartDistToBinXML.R")
#  doConversion("ascii")
#  doConversion("base64")
#  doConversion("raw")
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

# Function to do the encoding
doEncoding <- function(partVec, encoding) {
  if (encoding == "base64") {
    partVec = memCompress(partVec, "gzip")
    partVec = base64(partVec)
  } else if (encoding == "raw") {
    partVec = memCompress(partVec, "gzip")
    partVec = base64(partVec, TRUE, "raw")
  }
  return(partVec)
}

# Function to do the conversion
#doConversion <- function(boundaryCSV, encoding = "ascii") {

  # Input boundary csv file
  boundaryCSV = "input4000_boundary.csv"
  encoding = "ascii"

  # Set up compression algorithm (only gzip for now)
  compression = "none"
  if (encoding == "base64") {
    compression = "gzip"
  } else if (encoding == "raw") {
    compression = "gzip"
  }

  #------------------------------------------------------------------
  # Create XML
  # <?xml version='1.0' encoding='ISO-8859-1' ?>
  #<Ellip3D_input>
  #  <Meta>
  #   <title>  Boundary title </title>
  #  </Meta>
  #  <Boundary>
  #    <containerMin>  [0.0,  -1.0,  0.0] </containerMin>
  #    <containerMax>  [100.0, 1.0, 50.0] </containerMax>
  #    <boundary type="plane" id="x-">
  #      <direction> [-1.0, 0.0, 0.0] </direction>
  #      <position> [0.0, 0.0, 25.0] </position>
  #      <extraEdge>
  #        <direction> [-1.0, 0.0, 0.0] </direction>
  #        <position> [100.0, 0.0, 25.0] </position>
  #      </extraEdge>
  #    </boundary>
  #    <boundary type="plane" id="x+">
  #      <direction> [1.0, 0.0, 0.0] </direction>
  #      <position> [100.0, 0.0, 25.0] </position>
  #    </boundary>
  #    <boundary type="plane" id="y-">
  #      <direction> [0.0, -1.0, 0.0] </direction>
  #      <position> [50.0, -1.0, 25.0] </position>
  #    </boundary>
  #    <boundary type="plane" id="y+">
  #      <direction> [0.0, 1.0, 0.0] </direction>
  #      <position> [50.0, 1.0, 25.0] </position>
  #    </boundary>
  #    <boundary type="plane" id="z-">
  #      <direction> [0.0, 0.0, -1.0] </direction>
  #      <position> [50.0, 0.0, 0.0] </position>
  #    </boundary>
  #    <boundary type="plane" id="z+">
  #      <direction> [0.0, 0.0, 1.0] </direction>
  #      <position> [50.0, 0.0, 50.0] </position>
  #    </boundary>
  #  </Boundary>
  #</Ellip3D_input>
  #------------------------------------------------------------------

  # Read the first line (contains the boundary type and number)
  # bc_type = 0 => rigid boundary
  # bc_type = 1 => flexible boundary
  df_BC = read.csv(boundaryCSV, header = FALSE, nrows = 1, sep="", 
                   blank.lines.skip = TRUE)
  bc_type = df_BC$V1
  numBoundary = df_BC$V2
  print(numBoundary)

  boundary_geom = array()
  boundary_id = array()
  boundary_limit = array()
  boundary_area = array()
  boundary_order = list()
  boundary_radius = list()
  boundary_side = list()
  boundary_direction = list()
  boundary_position = list()
  lineCount = 1;
  for (boundary in 1:numBoundary) {
    # Read the boundary geometry
    # bc_geom = 1 => plane boundary
    # bc_geom = 2 => cylinder boundary
    df_BC = read.csv(boundaryCSV, header = FALSE, nrows = 1, skip = lineCount, sep="", 
                     blank.lines.skip = TRUE)
    lineCount = lineCount + 1
    boundary_geom[boundary] = df_BC$V1
    print(boundary_geom)

    # Read the boundary id, num limits, and area
    df_BC = read.csv(boundaryCSV, header = FALSE, nrows = 1, skip = lineCount, sep="", 
                     blank.lines.skip = TRUE)
    lineCount = lineCount + 1
    boundary_id[boundary] = df_BC$V1
    boundary_limit[boundary] = df_BC$V2
    boundary_area[boundary] = df_BC$V3
    print(paste(boundary_id, boundary_limit[boundary], boundary_area))

    # Read the coordinates of the limits
    df_BC = read.csv(boundaryCSV, header = FALSE, nrows = boundary_limit[boundary], 
                     skip = lineCount, sep="", blank.lines.skip = TRUE)
    lineCount = lineCount + boundary_limit[boundary] + 1
    print(df_BC)
    boundary_order[[boundary]] = df_BC$V1
    boundary_direction[[boundary]] = cbind(df_BC$V2, df_BC$V3, df_BC$V4)
    boundary_position[[boundary]] = cbind(df_BC$V5, df_BC$V6, df_BC$V7)
    boundary_radius[[boundary]] = df_BC$V8
    boundary_side[[boundary]] = df_BC$V9
    print(boundary_order)
    print(boundary_direction)
    print(boundary_position)
    print(boundary_radius)
    print(boundary_side)
  }

  # Compute the min and max of the boundary box (assuming axis-aligned)
  minX = 1.0e20
  minY = 1.0e20
  minZ = 1.0e20
  maxX = -1.0e20
  maxY = -1.0e20
  maxZ = -1.0e20
  for (boundary in 1:numBoundary) {
    minX = min(minX, min(boundary_position[[boundary]][,1]))
    minY = min(minY, min(boundary_position[[boundary]][,2]))
    minZ = min(minZ, min(boundary_position[[boundary]][,3]))
    maxX = max(maxX, max(boundary_position[[boundary]][,1]))
    maxY = max(maxY, max(boundary_position[[boundary]][,2]))
    maxZ = max(maxZ, max(boundary_position[[boundary]][,3]))
  }
  print(paste(minX, minY, minZ))
  print(paste(maxX, maxY, maxZ))
  

  # Create XML
  xml <- xmlOutputDOM(tag = "Ellip3D_input")

  xml$addTag("Meta", close = FALSE)
  xml$addTag("title", "Boundary file generated with convertBoundaryToXML.R" )
  xml$closeTag() # Close the Meta tag

  xml$addTag("Boundary", close = FALSE)

  #    <containerMin>  [0.0,  -1.0,  0.0] </containerMin>
  #    <containerMax>  [100.0, 1.0, 50.0] </containerMax>
  xml$addTag("containerMin", paste("[", minX, ",", minY, ",", minZ, "]"))
  xml$addTag("containerMax", paste("[", maxX, ",", maxY, ",", maxZ, "]"))

  #    <boundary type="plane" id="x-">
  #      <direction> [-1.0, 0.0, 0.0] </direction>
  #      <position> [0.0, 0.0, 25.0] </position>
  #      <extraEdge>
  #        <direction> [-1.0, 0.0, 0.0] </direction>
  #        <position> [100.0, 0.0, 25.0] </position>
  #      </extraEdge>
  #    </boundary>
  for (boundary in 1:numBoundary) {
    if (boundary_geom[boundary] == 1) {
      geometry = "plane"
    } else {
      geometry = "cylinder"
    }
    dir = boundary_direction[[boundary]][1,]
    if (all(dir == c(-1, 0, 0))) {
      direction = "x-"
    } else if (all(dir == c(1, 0, 0))) {
      direction = "x+"
    } else if (all(dir == c(0, -1, 0))) {
      direction = "y-"
    } else if (all(dir == c(0, 1, 0))) {
      direction = "y+"
    } else if (all(dir == c(0, 0, -1))) {
      direction = "z-"
    } else if (all(dir == c(0, 0, 1))) {
      direction = "z+"
    }
    pos = boundary_position[[boundary]][1,]
    print(paste(boundary, 1, geometry, direction, "dir:",
                dir[1], dir[2], dir[3], "pos:", pos[1], pos[2], pos[3]))
    xml$addTag("boundary", attrs = c(type = geometry, id = direction), close = FALSE)
    xml$addTag("direction", paste("[", dir[1], ",", dir[2], ",", dir[3], "]"))
    xml$addTag("position", paste("[", pos[1], ",", pos[2], ",", pos[3], "]"))
    if (bc_type == 0) {
      xml$addTag("initial_velocity", paste("[", 0, ",", 0, ",", 0, "]"))
    }

    if (boundary_limit[boundary] > 1) {
      for (jj in 2:boundary_limit[boundary]) {
        xml$addTag("extraEdge", close = FALSE)
        dir = boundary_direction[[boundary]][jj,]
        pos = boundary_position[[boundary]][jj,]
        print(paste(boundary, jj, "dir:",
                    dir[1], dir[2], dir[3], "pos:", pos[1], pos[2], pos[3]))
        xml$addTag("direction", paste("[", dir[1], ",", dir[2], ",", dir[3], "]"))
        xml$addTag("position", paste("[", pos[1], ",", pos[2], ",", pos[3], "]"))
        xml$closeTag() # Close the extraEdge tag
      }
    }
    xml$closeTag() # Close the boundary tag
  }

  xml$closeTag() # Close the Boundary tag

  #------------------------------------------------------------------
  # Save the XML file
  #------------------------------------------------------------------
  boundaryXML = ""
  if (encoding == "ascii") {
    boundaryXML = paste0(sub(".csv", "", boundaryCSV), ".ascii.xml")
  } else if (encoding == "base64") {
    boundaryXML = paste0(sub(".csv", "", boundaryCSV), ".base64.xml")
  } else {
    boundaryXML = paste0(sub(".csv", "", boundaryCSV), ".raw.xml")
  }
  cat(xmlFormat(xml), sep = "\n", file = boundaryXML)
# Function to do the conversion
doConversion <- function(boundaryCSV, encoding = "ascii") {

  startID = 1
  partShapes = unique(df_part$type)
  for (partShape in partShapes) {
    partShapeName = "ellipsoid"
    if (partShape == 1) {
      partShapeName = "sphere"
    } else if (partShape == 2) {
      partShapeName = "cube"
    } else if (partShape == 3) {
      partShapeName = "polyellipsoid"
    } else {
      partShapeName = "ellipsoid"
    }
    df = df_part[which(df_part$type == partShape),]

    # Add Id
    # Create a sequence of unique particle ids
    part_ids = seq(from = startID, to = (startID + nrow(df) - 1), by = 1)
    #partID = paste(df$id, collapse=" ")
    partID = paste(part_ids, collapse=" ")
    partID = doEncoding(partID, encoding)
    xml$addTag("id", partID, attrs = c(unit = "none", numComponents = "1"))

    # Add Type
    partType = paste(df$type, collapse=" ")
    partType = doEncoding(partType, encoding)
    xml$addTag("type", partType, attrs = c(unit = "none", numComponents = "1"))

    # Add radius
    radius = cbind(df$radius_a, df$radius_b, df$radius_c)
    partVec = paste(t(radius), collapse=" ")
    partVec = doEncoding(partVec, encoding)
    xml$addTag("radii", partVec, attrs = c(unit = "m", numComponents = "3"))

    # Add axle_a
    axle_a = cbind(df$axle_a_x, df$axle_a_y, df$axle_a_z)
    partVec = paste(t(axle_a), collapse=" ")
    partVec = doEncoding(partVec, encoding)
    xml$addTag("axle_a", partVec, attrs = c(unit = "rad", numComponents = "3"))

    # Add axle_b
    axle_b = cbind(df$axle_b_x, df$axle_b_y, df$axle_b_z)
    partVec = paste(t(axle_b), collapse=" ")
    partVec = doEncoding(partVec, encoding)
    xml$addTag("axle_b", partVec, attrs = c(unit = "rad", numComponents = "3"))

    # Add axle_c
    axle_c = cbind(df$axle_c_x, df$axle_c_y, df$axle_c_z)
    partVec = paste(t(axle_c), collapse=" ")
    partVec = doEncoding(partVec, encoding)
    xml$addTag("axle_c", partVec, attrs = c(unit = "rad", numComponents = "3"))

    # Add position
    position = cbind(df$position_x, df$position_y, df$position_z)
    partVec = paste(t(position), collapse=" ")
    partVec = doEncoding(partVec, encoding)
    xml$addTag("position", partVec, attrs = c(unit = "m", numComponents = "3"))

    # Add velocity
    velocity = cbind(df$velocity_x, df$velocity_y, df$velocity_z)
    partVec = paste(t(velocity), collapse=" ")
    partVec = doEncoding(partVec, encoding)
    xml$addTag("velocity", partVec, attrs = c(unit = "m/s", numComponents = "3"))

    # Add omega
    omega = cbind(df$omga_x, df$omga_y, df$omga_z)
    partVec = paste(t(omega), collapse=" ")
    partVec = doEncoding(partVec, encoding)
    xml$addTag("omega", partVec, attrs = c(unit = "rad/s", numComponents = "3"))

    # Add force
    force = cbind(df$force_x, df$force_y, df$force_z)
    partVec = paste(t(force), collapse=" ")
    partVec = doEncoding(partVec, encoding)
    xml$addTag("force", partVec, attrs = c(unit = "N", numComponents = "3"))

    # Add moment
    moment = cbind(df$moment_x, df$moment_y, df$moment_z)
    partVec = paste(t(moment), collapse=" ")
    partVec = doEncoding(partVec, encoding)
    xml$addTag("moment", partVec, attrs = c(unit = "Nm", numComponents = "3"))

    xml$closeTag() # Close the Particles tage

    startID = startID + nrow(df)
  }

  #------------------------------------------------------------------
  # Save the XML file
  #------------------------------------------------------------------
  boundaryXML = ""
  if (encoding == "ascii") {
    boundaryXML = paste0(sub(".csv", "", boundaryCSV), ".ascii.xml")
  } else if (encoding == "base64") {
    boundaryXML = paste0(sub(".csv", "", boundaryCSV), ".base64.xml")
  } else {
    boundaryXML = paste0(sub(".csv", "", boundaryCSV), ".raw.xml")
  }
  cat(xmlFormat(xml), sep = "\n", file = boundaryXML)
}

#doConversion("input4000_boundary.csv", "ascii")


