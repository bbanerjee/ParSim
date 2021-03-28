---
title:  "Reading VTK particles in Javascript"
subheadline: "Biswajit Banerjee"
description: "Parsing simulation output data"
date:  2017-02-21 09:30:00
categories:
    - Javascript
    - Typescript
    - Vue
    - Vuex
    - XML
excerpt_separator: <!--more-->
toc: true
toc_label: "Contents"
toc_sticky: true
---

In our [previous article](http://www.parresianz.com/c++/xml/vtk/vtk-particle-output/)
we discussed the process of calling XML functions in C++ to write data produced by
particle simulations in XML format files.  These output files can then
be viewed in powerful tools such as VisIt or Paraview.
<!--more-->

Modern browsers have become powerful in recent years and can now render
3D graphics quite efficiently.  The next few blog posts will discuss our
experience in developing a basic user interface with standard web tools
to visualize simulation output data.  In this article, I will explain
how the input process would work in such a tool.

##### The VTK XML particle data format #####
Output data from our simulations are written to disk every few timesteps.  At
each timestep, we get a data file like the one listed below.

{% highlight xml %}
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
  <UnstructuredGrid>
    <FieldData>
      <DataArray type="Float64" Name="TIME" NumberOfTuples="1" format="ascii" RangeMin="0" RangeMax="0">
        0.001
      </DataArray>
    </FieldData>
    <Piece NumberOfPoints="412" NumberOfCells="0">
      <PointData>
        <DataArray type="Float64" Name="Radius" NumberOfComponents="3" format="ascii" RangeMin="0.83645083538" RangeMax="1.7075128111">
          1.15 0.5 0.69 1.15 0.5 0.69
          1.15 0.5 0.69 0.575 0.5 0.345
          1.15 0.5 0.69 0.575 0.5 0.345
          0.575 0.5 0.345 1.15 0.5 0.69
          ...........
        </DataArray>
        <DataArray type="Float64" Name="Axis a" NumberOfComponents="3" format="ascii" RangeMin="1.923869419" RangeMax="3.8462926819">
          2.250054 1.570796 0.6792575 2.536687 1.570796 0.9658907
          0.8923825 1.570796 0.6784139 1.576513 1.570796 0.005716666
          1.523425 1.570796 0.0473713 2.090659 1.570796 2.62173
          2.128077 1.570796 0.5572808 0.0678937 1.570796 1.502903        
          ...........
        </DataArray>
        ................................................
      </PointData>
      <CellData>
      </CellData>
      <Points>
        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii" RangeMin="29.109286806" RangeMax="99.675090711">
          49.731739044 0 0.99367249012 50.900959015 0 5.3205289841
          30.811420441 0 3.8980329037 54.012958527 0 5.2524261475
          71.06136322 0 7.7233600616 65.397392273 0 6.3970627785
          42.516990662 0 11.963620186 48.560340881 0 2.5312030315
          ...........
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">
        </DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii" RangeMin="1e+299" RangeMax="-1e+299">
        </DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>

{% endhighlight %}

For now, let us consider only data that have been written in ASCII format.

##### Reading the particle XML data #####
We will read in this data into a code that provides a user interface to visualize
the particle data.

For now, the details of the user interface are not important,
but it is worth noting the following : <br>
* it is written in [Typescript](https://www.typescriptlang.org/) <br>
* it uses the [Vue](https://vuejs.org/) framework for user interactions and [Vuex](https://vuex.vuejs.org/en/) for data management <br>
* Typescript type definitions are from [npm](https://www.npmjs.com/) and [av-ts](https://github.com/HerringtonDarkholme/av-ts) <br>
* the code is transpiled to [ES6](http://es6-features.org/#Constants)
  using `tsc` and [babel](https://babeljs.io/), and <br>
* packaged using [webpack](https://webpack.github.io/) before it is run on a browser
{: .notice}

First let us examine the actual reader code:
{% highlight js %}
// Function: parseAndSaveVTKXML
// Input:   xmlDoc : XMLDocument (The current version of jquery.d.ts declares the return type of parseXML as "any")
// Store:   pointData : any {} (A key-value object containing arrays of particle data)
public parseAndSaveVTKXML(xmlDoc : any) {
    // Use jquery to read the xmlDoc
    let $xml = $(xmlDoc);
    // Check that the file is actual a VTK XML file
    let fileType = $xml.find("VTKFile").attr('type');
    if (fileType != "UnstructuredGrid") {
      console.log("Invalid file type" + fileType);
      return;
    }
    // Read time stamp after trimming the string
    let timeStr = $xml.find("FieldData").find("DataArray").text().trim();
    let time = parseFloat(time);
    // Read the number of points in the data set
    let numPtsStr = $xml.find("Piece").attr('NumberOfPoints').trim();
    let numPts = parseFloat(numPtsStr);
    // Read the state data one by one
    var pointData:any = {};
    $($xml.find("PointData").find("DataArray")).each(
      function () {
        let data = $(this);  // "this" is now pointing to the data XMLDocument
        let key = data.attr("Name");  // Use the name of the variable as the key
        let type = data.attr("type");
        let numComponents = data.attr("NumberOfComponents"); 
        if (numComponents) {  // More than one component
          let numComp = parseFloat(numComponents);
          let dataArray = convertTo1DArray(data);
          let arrayOfVec : any = [];
          while (dataArray.length) {
            arrayOfVec.push(dataArray.splice(0,numComp));  // Convert into 2D array        
          }
          pointData[key] = arrayOfVec; // Save the data                      
        } else { // Only one component
          pointData[key] = convertTo1DArray(data);
        }
      }
    );

    // Read the position data
    let points = $xml.find("Points").find("DataArray");
    let numComponents = points.attr("NumberOfComponents");
    let dataArray = convertTo1DArray(points);
    let numComp = parseFloat(numComponents);
    let arrayOfVec : any = [];
    while (dataArray.length) {
      arrayOfVec.push(dataArray.splice(0,numComp));                   
    }
    pointData["Position"] = arrayOfVec;
    // Save the data
    Store.commit('SET_VTK_PARTICLE_DATA', pointData); // Save to a Vuex store
  }

{% endhighlight %}

The `convertTo1DArray` method has the form

{% highlight js %}
private convertTo1DArray(data : string) : number[] {
  let dataArray =
    data.text().trim()                 // Remove leading/trailing white spaces
        .replace( /\s\s+/g, ' ' )      // Replace multiple spaces with one space
        .replace(/(\r\n|\n|\r)/gm," ") // Replace newline characters with space
        .split(" ")                    // Convert into 1D array
        .map(parseFloat);              // Convert into floats
  return dataArray;
}
{% endhighlight %}

Note that we have not saved the time in this simplified version. The JSON representation
of the saved data is of the form

{% highlight json %}
{
  "Position": [[1,2,3],[4,5,6],....],
  "Radius": [[0.1,0.2,0.3],[0.4,0.5,0.6],....],
  ....
}
{% endhighlight %}

##### Asynchronous IO and reading the file #####
The previous section assumes that a file has been read in and parsed using `parseXML`.
Let us now see how the file is actually read from disk.  We will assume that a
[File object](https://developer.mozilla.org/en-US/docs/Web/API/File) is available and
has been obtained using the `<input>` tag in the HTML docsument. The `FileReader` feature
of HTML5 makes it possible to read files easily from disk without using the `fs` library from
`node.js`.

The HTML input tag inside the appropriate `Vue` template has the form

{% highlight html %}
  <input type="file" v-on:change="readVTKXMLParticleFile">
{% endhighlight %}

The `readVTKXMLParticleFile` function is called after a file has been selected from the
list of files.

{% highlight js %}
public readVTKXMLParticleFile(event: any) {
  // Choose first file in selected list
  let file = event.target.files[0];
  // Read the file as text and store data
  this.readXMLFile(file, "vtk-xml");
}
{% endhighlight %}

The actual asynchronous read is defined in `readXMLFile`:
{% highlight js %}
import $ = require("jquery"); // Use jQuery for the parsing
.....
  public readXMLFile(file: any, type: string) {

    if (file) {
      let reader = new FileReader();
      reader.readAsText(file);
      // Use a closure to pass on the actual parser function
      reader.onload = (
        function(parserFunction: any) {
          return function(event: any) {
            let xmlData = reader.result;      // The text that has been read in
            let xmlDoc = $.parseXML(xmlData); // Use jQuery to do the parsing
            parserFunction(xmlDoc);
          };
        }
      )(this.parseAndSaveVTKXML); // This is where the actual parsing is done

    } else {
      console.log("Unable to load file ", file.name);
    }
  }
{% endhighlight %}

#### Remarks ####
Once the data have been stored, we can start the process of visualization.  Future articles will
discuss our experience with [vtk.js](https://kitware.github.io/vtk-js/) and [three.js](https://threejs.org/).  We will also discuss the potential and pitfalls of Typescript and Vue.

