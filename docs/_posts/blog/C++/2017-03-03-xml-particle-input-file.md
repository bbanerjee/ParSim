---
layout: posts
title:  "XML format for particle input files"
subheadline: "Biswajit Banerjee"
description: ""
date:  2017-03-03 10:30:00
categories:
    - R
    - XML
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

##### Introduction #####
Let us now change tack slightly and return to an issue I had talked about earlier in [Input files: reading XML]({{ site.baseurl }}/c++/xml/xml-input/).  Typical input files in research simulation codes cannot be easily deciphered.  But in some cases headers are included to make the reading process easier.

An example is a particle input file for a DEM code, called `particle_distribution.csv`, that looks
like the following:

~~~
 12
             id           type       radius_a       radius_b       radius_c     position_x     position_y     position_z       axle_a_x       axle_a_y       axle_a_z       axle_b_x       axle_b_y       axle_b_z       axle_c_x       axle_c_y       axle_c_z     velocity_x     velocity_y     velocity_z         omga_x         omga_y         omga_z        force_x        force_y        force_z       moment_x       moment_y       moment_z
              1              0   1.150000e+00   5.000000e-01   6.900000e-01   4.973174e+01   0.000000e+00   9.936725e-01   2.250054e+00   1.570796e+00   6.792575e-01   1.570796e+00   4.619360e-07   1.570796e+00   2.462335e+00   1.570796e+00   2.250054e+00   2.008031e-02   0.000000e+00   5.505538e-02   0.000000e+00   2.781917e-02   0.000000e+00  -4.916840e+05   0.000000e+00  -1.269380e+05   1.688596e-01  -1.187120e+05  -3.092984e-01
            187              0   1.150000e+00   5.000000e-01   6.900000e-01   5.090096e+01   0.000000e+00   5.320529e+00   2.536687e+00   1.570796e+00   9.658907e-01   1.570797e+00   1.803533e-06   1.570798e+00   2.175704e+00   1.570796e+00   2.536685e+00  -3.014865e-01   0.000000e+00  -1.031938e-01   0.000000e+00  -7.066198e-02   0.000000e+00  -4.413702e+03   0.000000e+00  -1.080548e+04  -1.174683e-01  -3.673358e+04  -5.749536e-02
            ........................
              5
   1.000000e+00   1.400000e+00
   8.000000e-01   1.300000e+00
   6.000000e-01   1.200000e+00
   3.000000e-01   1.150000e+00
   1.000000e-01   1.000000e+00

   8.000000e-01   6.000000e-01
~~~

There are 12 particles in this file and at the end there are some sieve sizes and passing
percentages, and also a radius ratio.  We would like to automatically convert this file into
a form that's more easily read and which also contains information about units. `R` provides
an excellent set of tools for that type of automation.

We would like to convert this data into a form that's easier to understand, e.g., 

{% highlight xml %}
  <?xml version='1.0'?>
  <Ellip3D_input>
    <Particles number="100" type="ellipsoid" compression="gzip" encoding="base64">
      <id unit="none" numComponents="1"> 1 2 3 ... 100 </id>
      <radii unit="m" numComponents="3" > a1 b1 c1 a2 b2 c2 ... </radii>
      <position unit="m" numComponents="3" > x1 y1 z1 x2 y2 z2 ... </position>
      <axle_a_angle unit="rad" numComponents="3" > theta_x1 theta_y1 theta_z1 ... </axle_a_angle>
      <axle_b_angle unit="rad" numComponents="3" > theta_x1 theta_y1 theta_z1 ... </axle_b_angle>
      <axle_c_angle unit="rad" numComponents="3"> theta_x1 theta_y1 theta_z1 ... </axle_c_angle>
      <velocity unit="m/s" numComponents="3"> x1 y1 z1 x2 y2 z2 ... </velocity>
      <omega unit="rad/s" numComponents="3"> x1 y1 z1 x2 y2 z2 ... </omega>
      <force unit="N" numComponents="3"> x1 y1 z1 x2 y2 z2 ... </force>
      <moment unit="Nm" numComponents="3"> x1 y1 z1 x2 y2 z2 ... </moment>
    </Particles>
    <Sieves number="5" compression="none" encoding="ascii">
      <percent_passing> 0.1 0.2 0.3 ... <percent_passing>
      <size unit="mm"> 0.1 0.2 0.3 ... <size>
      <sieve_ratio>
        <ratio_ba>0.8</ratio_ba>
        <ratio_ca>0.6</ratio_ca>
      </sieve_ratio>
    </Sieves>
  </Ellip3D_input>
{% endhighlight %}

##### R script converter from CSV to XML ####
Here's a script called `convertParDistToBinXML.R` that does the conversion:

{% highlight r %}
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
#  doConversion("raw")
#--------------------------------------------------------------------------
# Set the working directory
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
doConversion <- function(partDistCSV, encoding = "ascii") {

  # Input particle distribution csv file
  #partDistCSV = "particle_distribution.csv"

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
  #  <Particles number="100" type="ellipsoid" compression="gzip" encoding="base64">
  #    <id unit="none" numComponents="1"> 1 2 3 ... 100 </id>
  #    <radii unit="m" numComponents="3" > a1 b1 c1 a2 b2 c2 ... </radii>
  #    <position unit="m" numComponents="3" > x1 y1 z1 x2 y2 z2 ... </position>
  #    <axle_a_angle unit="rad" numComponents="3" > theta_x1 theta_y1 theta_z1 ... </axle_a_angle>
  #    <axle_b_angle unit="rad" numComponents="3" > theta_x1 theta_y1 theta_z1 ... </axle_b_angle>
  #    <axle_c_angle unit="rad" numComponents="3"> theta_x1 theta_y1 theta_z1 ... </axle_c_angle>
  #    <velocity unit="m/s" numComponents="3"> x1 y1 z1 x2 y2 z2 ... </velocity>
  #    <omega unit="rad/s" numComponents="3"> x1 y1 z1 x2 y2 z2 ... </omega>
  #    <force unit="N" numComponents="3"> x1 y1 z1 x2 y2 z2 ... </force>
  #    <moment unit="Nm" numComponents="3"> x1 y1 z1 x2 y2 z2 ... </moment>
  #  </Particles>
  #  <Sieves number="5" compression="none" encoding="ascii">
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
    xml$addTag("Particles", 
               attrs = c(number = nrow(df), type = partTypeName, 
                         compression = compression, encoding = encoding), 
               close = FALSE)

    # Add Id
    partID = paste(df$id, collapse=" ")
    partID = doEncoding(partID, encoding)
    xml$addTag("id", partID, attrs = c(unit = "none", numComponents = "1"))

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
  xml$addTag("Sieves", 
             attrs = c(number = numSieve, compression = "none", 
                       encoding = "ascii"), 
             close = FALSE)
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
  } else if (encoding == "base64") {
    partDistXML = paste0(sub(".csv", "", partDistCSV), ".base64.xml")
  } else {
    partDistXML = paste0(sub(".csv", "", partDistCSV), ".raw.xml")
  }
  cat(xmlFormat(xml), sep = "\n", file = partDistXML)
}

doConversion("particle_distribution.short.csv", "ascii")
doConversion("particle_distribution.short.csv", "base64")
doConversion("particle_distribution.short.csv", "raw")
doConversion("particle_distribution.csv", "ascii")
doConversion("particle_distribution.csv", "base64")
doConversion("particle_distribution.csv", "raw")
{% endhighlight %}

You can run this script in `R` by first loading it using `source("convertPartDistToBinXML.R")`
and then running  `doConversion("file_name.csv", "ascii")` or `doConversion("file_name.csv", "base64")`.
{:.notice--info}

##### The output ASCII XML file #####
The ASCII XML file produced by this script is `particle_distribution.ascii.xml` and contains

{% highlight xml %}

<?xml version="1.0"?>
<Ellip3D_input>
  <Particles number="3" type="ellipsoid" compression="none" encoding="ascii">
    <id unit="none" numComponents="1">1 187 373</id>
    <radii unit="m" numComponents="3">1.15 0.5 0.69 1.15 0.5 0.69 1.15 0.5 0.69</radii>
    <axle_a unit="rad" numComponents="3">2.250054 1.570796 0.6792575 2.536687 1.570796 0.9658907 0.8923825 1.570796 0.6784139</axle_a>
    <axle_b unit="rad" numComponents="3">1.570796 4.61936e-07 1.570796 1.570797 1.803533e-06 1.570798 1.570799 2.865908e-06 1.570795</axle_b>
    <axle_c unit="rad" numComponents="3">2.462335 1.570796 2.250054 2.175704 1.570796 2.536685 2.463179 1.570796 0.8923828</axle_c>
    <position unit="m" numComponents="3">49.73174 0 0.9936725 50.90096 0 5.320529 30.81142 0 3.898033</position>
    <velocity unit="m/s" numComponents="3">0.02008031 0 0.05505538 -0.3014865 0 -0.1031938 -0.1538927 0 0.04640176</velocity>
    <omega unit="rad/s" numComponents="3">0 0.02781917 0 0 -0.07066198 0 0 -0.2856931 0</omega>
    <force unit="N" numComponents="3">-491684 0 -126938 -4413.702 0 -10805.48 22966.5 0 6613.044</force>
    <moment unit="Nm" numComponents="3">0.1688596 -118712 -0.3092984 -0.1174683 -36733.58 -0.05749536 0.0414009 -5030.464 -0.0190764</moment>
  </Particles>
  <Particles number="3" type="sphere" compression="none" encoding="ascii">
    <id unit="none" numComponents="1">373 373 559</id>
    <radii unit="m" numComponents="3">0.575 0.5 0.345 1.15 0.5 0.69 0.575 0.5 0.345</radii>
    <axle_a unit="rad" numComponents="3">1.576513 1.570796 0.005716666 1.523425 1.570796 0.0473713 2.090659 1.570796 2.62173</axle_a>
    <axle_b unit="rad" numComponents="3">1.570798 1.803533e-06 1.570798 1.570798 2.756634e-06 1.570799 1.570794 2.349616e-06 1.570795</axle_b>
    <axle_c unit="rad" numComponents="3">3.135876 1.570796 1.576513 3.094221 1.570796 1.523425 0.519865 1.570796 2.090661</axle_c>
    <position unit="m" numComponents="3">54.01296 0 5.252426 71.06136 0 7.72336 65.39739 0 6.397063</position>
    <velocity unit="m/s" numComponents="3">-0.3850849 0 -0.1090706 -0.173833 0 -0.03869174 -0.08703834 0 -0.02640878</velocity>
    <omega unit="rad/s" numComponents="3">0 0.1298346 0 0 0.06904687 0 0 0.271116 0</omega>
    <force unit="N" numComponents="3">120004 0 -76384.24 192576.1 0 108463.9 83898.65 0 115703.8</force>
    <moment unit="Nm" numComponents="3">-0.01118079 -6253.187 0.1382343 -0.2065214 6178.625 0.2630611 0.1774521 2937.92 -0.3412123</moment>
  </Particles>
  .....
  <Sieves number="5" compression="none" encoding="ascii">
    <percent_passing>1 0.8 0.6 0.3 0.1</percent_passing>
    <size unit="mm">1.4 1.3 1.2 1.15 1</size>
    <sieve_ratio>
      <ratio_ba>0.8</ratio_ba>
      <ratio_ca>0.6</ratio_ca>
    </sieve_ratio>
  </Sieves>
</Ellip3D_input>
{% endhighlight %}

##### The output compressed base64 XML file #####
The compressed base64 XML file produced by this script is `particle_distribution.base64.xml` and contains

{% highlight xml %}
<?xml version="1.0"?>
<Ellip3D_input>
  <Particles number="3" type="ellipsoid" compression="gzip" encoding="base64">
    <id unit="none" numComponents="1">eJwzVDC0MFcwNjcGAAg2Aa8=</id>
    <radii unit="m" numComponents="3">eJwz1DM0VTDQA2EzSwVD3DwAm34HcA==</radii>
    <axle_a unit="rad" numComponents="3">eJxVy8ENwEAIA8FWqABx5oxx/40lz+S30miRYBVvnKRKnqgcGRQDyZ5Zfc3DdemtNXrB/7j3tB+22RBv</axle_a>
    <axle_b unit="rad" numComponents="3">eJxNyckNACAIAMFWaEACIldB9t+Cmojxt5llVCdPg4HGKTYbOXDhjSNBoiJ7P42KhI5hmhT/1gVJPxJZ</axle_b>
    <axle_c unit="rad" numComponents="3">eJxNyckNADAIxMBWUgGCheXov7HwS34jGxIJdx4TltbkgYCqjIXVtvgXPbO5iHSreUulB97oC4BQD7g=</axle_c>
    <position unit="m" numComponents="3">eJwNyMERAEAEA8BWVGBCONJ/Y+e3syUfxpTB4BLfZFsfAb3LdiY6ZYRvROUdfbUgP13AC20=</position>
    <velocity unit="m/s" numComponents="3">eJwlitsNACAIxFZhAUmRh7j/YpKY3M+1RdnQuAmCkjnzloU6Fl05fI5NcT+38Xef30cFduoB3dgNXA==</velocity>
    <omega unit="rad/s" numComponents="3">eJwzUDDQMzAytzC0NDRXMABCXSDf3MDMzNDSAsY3sjA1szQ2VDAAAL6CCE8=</omega>
    <force unit="N" numComponents="3">eJwdycERACAIA7BVWACuYK2w/2J6fhPnpJoG8yzNanMyVxzUNzR2sK1qpNjPpNcgL0K4Cu8=</force>
    <moment unit="Nm" numComponents="3">eJwdy8kNwDAMA8FW1IAE0rr7byyxf4MFFsaayS1RcppHFObYsxOXZEeNi3q1u+XciOzY9JKfwQBWNOGwqDeBi674ADeiEak=</moment>
  </Particles>
  <Particles number="3" type="sphere" compression="gzip" encoding="base64">
    <id unit="none" numComponents="1">eJwzNjdWMAZiU1NLAAydAh4=</id>
    <radii unit="m" numComponents="3">eJwz0DM1N1Uw0ANhYxNTBUM9QxjXzBLEQJYFAL1ICD4=</radii>
    <axle_a unit="rad" numComponents="3">eJxVi8ENwEAMwla5CSJCjqDsv1jTV1Veli0y5FbWyQV4+iAAOXv3StalfvW6vAcGBq35GqOZrgeh0xAM</axle_a>
    <axle_b unit="rad" numComponents="3">eJwz1DM1NzC3tFAw1LMwMDY1Nk7VNTADcuCiUIaRnrmpmZmxCbK0JYxhApQ2NrE0MzRDljYFAFruEoY=</axle_b>
    <axle_c unit="rad" numComponents="3">eJxVysENADEIA8FWqABhiE3ov7GcdI8ov5Vmy1HcLYOzo+cPEWXlMSsTD2WtpIUTs8VL+c0h4QByxA96</axle_c>
    <position unit="m" numComponents="3">eJwVyEENADAIA0ArKGigHSX4N7btd7k+yOI6MhpsHjqmkC79Gwz15IZ2tK/8ldYFQ+cLCg==</position>
    <velocity unit="m/s" numComponents="3">eJwtickRACAIA1uxAZ1oEEL/jemgvz06BrUhy4bWMSYSAS8MinwZlOcMK1Zcpf2z3G7RAfzcDdI=</velocity>
    <omega unit="rad/s" numComponents="3">eJwzUDDQMzSytDA2MVMwAEE9AzNLAxMzC3Mo18jc0NAQKAcAnmMHhg==</omega>
    <force unit="N" numComponents="3">eJwVyEkRwAAIA0ArGCiTcAb/xtrucxkAymDPTqo8ynjRO84vCdWknyl18un/2It0vUKNCwM=</force>
    <moment unit="Nm" numComponents="3">eJwVysERACAIA7BVWECOtgi4/2LqO1nhAWCij63ilmPawqGhUrbCGbWJtEKPv/KUpSjgv+58ajxqP/xfCYK6DRMRKw==</moment>
  </Particles>
  .....
  <Sieves number="5" compression="none" encoding="ascii">
    <percent_passing>1 0.8 0.6 0.3 0.1</percent_passing>
    <size unit="mm">1.4 1.3 1.2 1.15 1</size>
    <sieve_ratio>
      <ratio_ba>0.8</ratio_ba>
      <ratio_ca>0.6</ratio_ca>
    </sieve_ratio>
  </Sieves>
</Ellip3D_input>
{% endhighlight %}

For a file containing 412 particles, the original CSV file was of size 177K, the ASCII XML version was 85K, and
the compressed base64 XML file was 39.  The compressed ASCII XML file was also smaller than the compressed base64 XML file.
{:.notice--info}

#### Remarks ####
In general it is better to use structured files for input data, provided the size of the data
is not too large.  If the input files are large, a compressed binary format
is better and should be used if possible.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

