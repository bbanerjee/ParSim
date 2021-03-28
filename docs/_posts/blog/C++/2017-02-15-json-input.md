---
layout: posts
title:  "Reading JSON in C++"
subheadline: "Biswajit Banerjee"
description: "JSON input files for research codes"
date:  2017-02-15 09:30:00
categories:
    - C++
    - JSON
excerpt_separator: <!--more-->
toc: true
toc_label: "Contents"
toc_sticky: true
---

In my previous post on [reading XML input files](http://www.parresianz.com/c++/xml/xml-input/), I discussed how input files can be made more human friendly with XML markup.  In some situations, it may be preferable to have [JSON format files](http://json.org/example.html) instead.
<!--more-->

JSON is particularly useful when the [Javascript](http://es6-features.org/#Constants) is used during reading and writing.  Though Javascript is not the ideal platform for computational engineering, it has some potential is used for user interface development.

*Note*: It remains to be seen how useful Javascript user interfaces are when used for engineering applications, though cross-platform tools such as [Electron](http://electron.atom.io/) has potential.  User interfaces designed with [wxWidgets](https://www.wxwidgets.org/) or [Qt](https://www.qt.io/) still dominate the landscape for engineering software.
{: .notice}

##### The JSON input file #####
In JSON format, our `Boundary` input file can be expressed as:
{% highlight json %}
  {"Boundary": {
    "containerMin":  "[0.0,  -1.0,  0.0]",
    "containerMax":  "[100.0, 1.0, 50.0]",
    "boundary": [
      {
       "type": "plane", "id": 1,
       "direction": "[-1.0, 0.0, 0.0]",
       "position": "[100.0, 0.0, 25.0]",
      },
      {
       "type": "plane", "id": 2,
       "direction": "[1.0, 0.0, 0.0]",
       "position": "[100.0, 0.0, 25.0]"
      }
    ]
  }}
{% endhighlight %}

Recall that the equivalent XML input file  has the form
{% highlight xml %}
  <?xml version='1.0' encoding='ISO-8859-1' ?>
  <Boundary>
    <!-- Container limits -->
    <containerMin>  [0.0,  -1.0,  0.0] </containerMin>
    <containerMax>  [100.0, 1.0, 50.0] </containerMax>
    <!-- Internal boundaries -->
    <boundary type="plane" id="1">
      <direction> [-1.0, 0.0, 0.0] </direction>
      <position> [100.0, 0.0, 25.0] </position>
    </boundary>
    <boundary type="plane" id="2">
      <direction> [1.0, 0.0, 0.0] </direction>
      <position> [100.0, 0.0, 25.0] </position>
    </boundary>
  </Boundary>
{% endhighlight %}

I prefer the XML version, but most new students will be familiar with
JSON and may prefer that format, particularly as that format appears to
be preferred in IoT applications.

##### Reading JSON #####
A large list of JSON readers for C++ can be found in the [JSON webdocs](http://www.json.org/).  A reasonably good header-only library is [JSON for Modern C++](https://github.com/nlohmann/json).  I have chosen to use that library because its constructs are similar to those of [ZenXml](http://zenxml.sourceforge.net/).

We add a new method, `readJSON`, to the class `BoundaryReader` and
implement it as follows.

{% highlight cpp %}
#include <BoundaryReader.h>
#include <json.hpp>
using json = nlohmann::json;
bool
BoundaryReader::readXML(const std::string& inputFileName) const 
{
  // Create an input filestream
  std::ifstream ifs(inputFileName);
  if (!ifs) {
    std::cout << "*ERROR** Could not read input file " << inputFileName << "\n";
    return false;
  }
  // Read the file stream into a string stream
  std::stringstream iss;
  iss << ifs.rdbuf();
  ifs.close();
  // Parse the input stream
  json docs;
  try {
    docs << iss;
  } catch (std::invalid_argument e) {
    std::cout << "*ERROR** Could not parse input file " << inputFileName
              << "\n";
    std::cout << "Please check for correctness using a linter.\n";
    return false;
  }
  // Read the boundary information
  json boundary_ps;
  try {
    boundary_ps = ps["Boundary"];
  } catch (std::exception e) {
    std::cout << "**ERROR** \"Boundary\" key not found. \n";
    return false;
  }
  // Read the container dimensions
  std::string vecStr;
  try {
    vecStr = boundary_ps["containerMin"].get<std::string>();
  } catch (std::exception e) {
    std::cout
      << "**ERROR** Container min. position not found in boundary geometry\n";
    std::cout << "  Add the containerMin: [x, y, z]  tag.";
    return false;
  }
  Vec boxMin = Vec::fromString(vecStr);
  try {
    vecStr = boundary_ps["containerMax"].get<std::string>();
  } catch (std::exception e) {
    .....
    return false;
  }
  Vec boxMax = Vec::fromString(vecStr);
  
  try {
    auto bound_ps = boundary_ps["boundary"];
    for (auto object : bound_ps) {
      std::string boundaryType = object["type"].get<std::string>();
      BoundaryId id = object["id"].get<BoundaryId>();
      vecStr = ps["direction"].get<std::string>();
      Vec direc = Vec::fromString(vecStr);
      vecStr = ps["position"].get<std::string>();
      Vec point = Vec::fromString(vecStr);
    }
  } catch (std::exception e) {
    std::cout << "**ERROR** Boundaries not found in boundary geometry\n";
    std::cout << "  Add the boundary key: value array tags etc.";
    return false;
  }
  return true;
}
{% endhighlight %}

I'm not really a fan of `try-catch` blocks, but that's the easiest way of dealing cleanly with read/parse failures when using `JSON for Modern C++`.  For research codes, this approach will suffice.  An alternative is to use a check of the form
{% highlight cpp %}
{
  auto ps_iter = ps.find("key_name");
  if (ps_iter == ps.end()) {
    // Error
  }
}
{% endhighlight %}

#### Remarks ####
Sometimes input/output files for engineering simulations can be quite large.  In those situations I would suggest XML with binary data. In future post I'll show you how to create a VTK XML format output file with binary data.


If you have questions/comments/corrections, please contact banerjee at parresianz dot com dot json (without the dot json).


