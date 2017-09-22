---
layout: docs
title:  "Reading XML files containing gzipped data in C++"
subheadline: "Biswajit Banerjee"
description: "How to read particle input files created with R in XML format"
date:  2017-04-04 10:30:00
categories:
    - C++
    - XML
    - gzip
image:
    credit: Parresia Research Limited
    header: "HummerLargeSim-WithLogo.png"
---

- Contents
{:toc}
{:.notice--content}

#### Introduction ####
We saw how to create an XML file containing compressed particle data in
the article ["XML format for particle input files"]({{site.url }}/r/xml/xml-particle-input-file/).
Let us now explore how to read in that data in our C++ particle simulation code.

#### Recap ####
Recall that the compressed base64 XML file contains data of the form shown below.  We would like
to convert these data back into numerical values that can be used by the simulation.

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

#### The header file ####
We use a `ParticleFileReader` object to read the file.  The declaration of the object is
listed below. Particle data are stored in an array of pointers to `Particle`
objects, called `ParticlePArray`.

{% highlight cpp %}
#include <zenxml/xml.h>
#include <string>
#include <vector>
class ParticleFileReader
{
public:
  ....
  void read(const std::string& fileName,
            ParticlePArray& particles) const;
private:
  template<typename T>
  bool readParticleValues(zen::XmlIn& ps,
                          const std::string& name,
                          const std::string& particleType,
                          std::vector<T>& output) const;
  template<typename T>
  bool decodeAndUncompress(const std::string &inputStr,
                           const int& numComponents,
                           std::vector<T>& output) const;
  template<typename T>
  T convert(const std::string& str) const;
  ....
};
{% endhighlight %}

The main workhorse methods in this class are the templated functions
`readParticleValues`, `decodeAndUncompress`, and `convert`.  Templates
are used because similar logic is used for different variable types.

#### The implementation ####
Let us now look at the implementations of these functions.  We will ignore
any checks that are necessary to make sure that the XML file is readable
and contains the right data.

##### The `read` function #####
The point of entry is the `read` function:

{% highlight cpp %}
bool
ParticleFileReader::read(const std::string &inputFileName,
                         ParticlePArray &particles) const
{
  // Read the input file
  zen::XmlDoc docs = zen::load(inputFileName);
  // Load the docsument into input proxy for easier element access
  zen::XmlIn ps(docs);
  // Loop through the particle types in the input file
  for (auto particle_ps = ps["Particles"]; particle_ps; particle_ps.next()) {
    // Get the attributes of the particles
    std::size_t numParticles = 0;
    std::string particleType = "sphere", compression = "none", encoding = "none";
    particle_ps.attribute("number", numParticles);
    particle_ps.attribute("type", particleType);
    particle_ps.attribute("compression", compression);
    particle_ps.attribute("encoding", encoding);
    // Assume that the input file is encoded and compressed
    if (encoding == "base64" && compression == "gzip") {
      // Get the particle ids
      std::vector<size_t> particleIDs;
      bool success = readParticleValues<size_t>(particle_ps,
                                                "id", particleType,
                                                particleIDs);
      // Get the particle radii
      std::vector<Vec> particleRadii;
      success = readParticleValues<Vec>(particle_ps,
                                        "radii", particleType,
                                        particleRadii);
      ............
      ............
      // Create the Particle array
      for (std::size_t ii = 0; ii < numParticles; ++ii) {
        ParticleP pt = std::make_shared<Particle>(
            particleIDs[ii], particleType, particleRadii[ii], ....);
        particles.push_back(pt);
      }
    } // end if (encoding == "base64" && compression == "gzip") 
  } // end for
  .....
}
{% endhighlight %}

The particle data associated with each tag is an array containing either 1 or 3
components.  We use the `explicitly instantiated` templated function
`readParticleValues<T>` to read in the data into arrays.

##### The `readParticleValues` templated function #####
Let us now look at the `readParticleValues` function that does the extraction
and conversion of the compressed and encoded data.

{% highlight cpp %}
template <typename T>
bool
ParticleFileReader::readParticleValues(zen::XmlIn &ps,
                                       const std::string &name,
                                       const std::string &particleType,
                                       std::vector<T> &output) const
{
  // Get the particle values
  int numComp = 1;
  std::string particleDataStr;
  auto prop_ps = ps[name];
  prop_ps.attribute("numComponents", numComp);
  prop_ps(particleDataStr);
  // Do the decoding and inflation of the compressed data
  bool success = decodeAndUncompress<T>(particleDataStr, numComp, output);
  return success;
}
{% endhighlight %}

The function just extracts the encoded data from the XML file and the number
of components in the data (1 or 3).  It then passes these on to the actual
decode and uncompress code.

##### The `decodeAndUncompress` templated function #####
This is where the main work is done.  For decoding the data into
binary form, we use the [cppcodec](https://github.com/tplgy/cppcodec)
library.  For decompression we use [ZLib](http://zlib.net/zlib_how.html).
To make sure that the `cppcodec` library is available in the repository
where our code is stored, we add it as a submodule using

~~~~
git submodule add git://github.com/tplgy/cppcodec.git cppcodec
~~~~

For the `Zlib` library to be available to our `cmake` build system, we add
the following to our `CMakeLists.txt` file:

~~~~
#-------------------------------------------------------
# Add requirements for Zlib compression library 
#-------------------------------------------------------
find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
  message(STATUS "Zlib compression library found")
  include_directories(${ZLIB_INCLUDE_DIRS})
else()
  message(STATUS "Zlib compression library not found")
  set(ZLIB_DIR "")
  set(ZLIB_LIBRARIES "")
  set(ZLIB_INCLUDE_DIRS "")
endif()
~~~~

The code for the `decodeAndUncompress` function is listed below.

{% highlight cpp %}
#include <cppcodec/cppcodec/base64_default_rfc4648.hpp>
#include "zlib.h"
template <typename T>
bool
ParticleFileReader::decodeAndUncompress(const std::string &inputStr,
                                        const int &numComponents,
                                        std::vector<T> &output) const
{
  // Decode from base64
  std::vector<std::uint8_t> decoded = base64::decode(inputStr);
  // Uncompress from gzip
  std::vector<std::uint8_t> uncompressed;
  z_stream stream;
  // Allocate inflate state
  stream.zalloc = Z_NULL;
  stream.zfree = Z_NULL;
  stream.opaque = Z_NULL;
  stream.avail_in = 0;
  stream.next_in = Z_NULL;
  int err = inflateInit(&stream);
  if (err != Z_OK) {
    std::cerr << "inflateInit" << " error: " << err << std::endl;
    return false;
  }
  // Uncompress until stream ends
  stream.avail_in = decoded.size();
  stream.next_in = &decoded[0];
  do {
    do {
      std::vector<std::uint8_t> out(decoded.size());
      stream.avail_out = out.size();
      stream.next_out = &out[0];
      err = inflate(&stream, Z_SYNC_FLUSH);
      uncompressed.insert(std::end(uncompressed), std::begin(out), std::end(out));
    } while (stream.avail_out == 0);
  } while (err != Z_STREAM_END);
  // Clean up and exit
  if (inflateEnd(&stream) != Z_OK) {
    std::cerr << "inflateEnd" << " error: " << err << std::endl;
    return false;
  }
  // Split the uncompressed string into a vector of tokens
  // (Assume that data are space separated)
  // (See: https://stackoverflow.com/questions/236129/split-a-string-in-c)
  std::istringstream iss(std::string(uncompressed.begin(), uncompressed.end()));
  std::vector<std::string> outputStr = {std::istream_iterator<std::string>{iss},
                                        std::istream_iterator<std::string>{}};
  // Convert the strings into the right type
  for (auto iter = outputStr.begin(); iter != outputStr.end(); iter += numComponents) {
    std::string str = *iter;
    for (int ii = 1; ii < numComponents; ii++) {
      // For more than one component, join into string with space separator
      str += " ";
      str += *(iter + ii);
    }
    output.push_back(convert<T>(str));
  }
  return true;
}
{% endhighlight %}

Here the main complication arises during the inflation of the compressed data.  We
don't know the size of the output buffer beforehand and have to read the buffer repeatedly
until the entire input buffer has been inflated.  After each chunk has been read
into the `out` vector, we insert the data into `uncompressed` and continue the process.

After the entire stream has be uncompressed, we convert the string into the correct
size type using the `convert<T>` function.  Notice that this function is `implicitly
instantiated` using `output.push_back(convert<T>(str))`.  Template specialization
is needed at this stage to make sure the right work is work during the conversion
of each type.  To see why this is not always a good idea, see the article 
[Why Not Specialize Function Templates?](http://www.gotw.ca/publications/mill17.htm).
Care is needed to make sure that we don't try to explictly instantiate `convert<T>`
elsewhere, and modern compilers will probably throw an error if that is attempted.

##### The `convert<T>` template specializations #####
We will define two specializations here;  the first function deals with
properties such as particle ID while the second deals with vector properties
such as position and force.

{% highlight cpp %}
template <>
size_t
ParticleFileReader::convert<size_t>(const std::string &str) const
{
  return std::stoul(str);
}
template <>
Vec
ParticleFileReader::convert<Vec>(const std::string &str) const
{
  std::istringstream iss(std::string(str.begin(), str.end()));
  std::vector<std::string> split = {std::istream_iterator<std::string>{iss},
                                    std::istream_iterator<std::string>{}};
  return Vec(std::stod(split[0]), std::stod(split[1]), std::stod(split[2]));
}
{% endhighlight %}

That completes the implementation.  To see a version of this approach in action,
look at [ParticleFileReader.cpp](https://github.com/bbanerjee/ParSim/blob/master/ThirdParty/paraEllip3d_DEM_PD/src/InputOutput/ParticleFileReader.cpp).

#### Remarks ####
We can see that the process or decoding and unzipping the data in the XML file is
quite straightforward.  But it takes a bit more effort than reading a formatted text file.
However, if our data include millions of particles, and these particles have to be
broadcast to several nodes of a multiprocessor system, compression can not only save us a
lot of communication time during simulations but also disk space.

In the next article, we will explore some more aspects of our particle simulation code.

If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot zen (without the dot zen).


<a class="twitter-share-button" href="https://twitter.com/intent/tweet" data-via="parresianz"> Tweet</a>
<script src="//platform.linkedin.com/in.js" type="text/javascript">
  lang: en_US
</script>
<script type="IN/Share" data-counter="right"></script>

<script src="https://d3js.org/d3.v4.min.js"></script>
<script src="{{ site.url }}/assets/js/yieldsurface.js"></script>
<script>
  d3.json("{{ site.url }}/assets/json/yieldSurfData.json", drawYieldSurface);
</script>

