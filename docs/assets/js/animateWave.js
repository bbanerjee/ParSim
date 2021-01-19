      "use strict";

      function drawCFLWaveAnimation() {

        // Create svg
        let svgsize = {x: 500, y: 300};
        let margin = {top: 50, right: 50, bottom: 50, left: 50};
        let width = svgsize.x - margin.left - margin.right;
        let height = svgsize.y - margin.top - margin.bottom;

        let svg = d3.select(".cfl-wave-animation")
                    .append("svg")
                      .attr("width", svgsize.x)
                      .attr("height", svgsize.y);
        let grid = svg.append("g")
                      .attr("class", "grid")
                      .attr("transform", 
                            "translate(" + margin.left + "," + margin.top + ")");

        // Create scaling functions
        let xmin = -6.0;
        let xmax =  6.0;
        let ymin =  0.0;
        let ymax =  5.0;
        let xscale = d3.scaleLinear().domain([xmin, xmax]).rangeRound([0, width]);
        let yscale = d3.scaleLinear().domain([ymin, ymax]).rangeRound([height, 0]);

        // Create polyline generation function
        let line = d3.line()
                     .x(function(d) { return xscale(d[0]); })
                     .y(function(d) { return yscale(d[1]); });

        // Create bottom axis
        grid.append("g")
          .attr("transform", "translate(0," + height + ")")
          .call(d3.axisBottom(xscale))
          .append("text")
          .attr("fill", "#000")
          .attr("transform",
                "translate(" + (width/2) + " ," + 
                             (margin.top - 20) + ")")
           .attr("text-anchor", "end")
           .text("Position (x)"); 

        // Create left axis
        grid.append("g")
          .call(d3.axisLeft(yscale))
          .append("text")
          .attr("fill", "#000")
          .attr("transform",
                "translate(" + (-margin.left) + " ," +
                                    (height/2) + ") rotate(-90)")
          .attr("y", 6)
          .attr("dy", "0.71em")
          .attr("text-anchor", "end")
          .text("Time (t)"); 

        // Plot horizontal lines
        let lines = [];
        for (let yCoord = ymin; yCoord < ymax+1; yCoord++) {
          lines[yCoord] = grid.append("g")
            .append("path")
            .attr("class", "line-y"+yCoord)
            .datum([[xmin, yCoord],[xmax, yCoord]])
            .attr("fill", "none")    
            .attr("stroke", "gray")
            .attr("stroke-width", 0.5)
            .attr("opacity", 0)
            .attr("d", line);
        }

        // Plot grid points
        let gridPoints = [];
        for (let yCoord = ymin; yCoord < ymax+1; yCoord++) {
          let gridData = [];
          let index = 0;
          for (let xCoord = xmin; xCoord < xmax+1; xCoord++) {
            gridData[index] = [xCoord, yCoord];
            index++;
          }
          gridPoints[yCoord] = grid.append("g")
            .attr("class", "grid-points"+yCoord)
            .attr("opacity", 0)
            .selectAll("grid-point-y"+yCoord)
            .data(gridData)
            .enter()
            .append("circle")
            .attr("r", 2)
            .attr("cx", function(d) {return xscale(d[0]);})
            .attr("cy", function(d) {return yscale(d[1]);})
            .attr("fill", "blue")
            .attr("stroke", "black");
        }
                                   
        // Plot wave front
        let waveFront = [];
        for (let yCoord = ymin; yCoord < ymax+1; yCoord++) {
          let waveFrontData = [];
          waveFrontData = [[(- yCoord), ymin],
                           [( yCoord), ymin],
                           [0, yCoord]]; 
          waveFront[yCoord] = grid.append("g")
            .append("path")
            .attr("class", "wave-y"+yCoord)
            .datum(waveFrontData)
            .attr("fill", "pink")    
            .attr("stroke", "gray")
            .attr("stroke-width", 0.5)
            .attr("opacity", 0)
            .attr("d", line);
        }
                

        function drawAnimation() {

          let yCoord = ymin;
          function drawWave() {
            d3.active(this)
              .attr("opacity", 1);
            grid.selectAll(".grid-points"+yCoord)
              .attr("opacity", 0.8);
            grid.select(".wave-y"+yCoord)
              .attr("opacity", 0.6);
            yCoord++;
            if (yCoord > ymax) {
              yCoord = ymin;
              resetWave();
            }
            grid.select(".line-y"+yCoord)
              .transition()
              .delay(1000)
              .duration(1000)
              .on("start", drawWave);
          }

          function resetWave() {
            for (let yy = ymin; yy < ymax+1; yy++) {
              grid.select(".line-y"+yy)
                .attr("opacity", 0.1);
              grid.selectAll(".grid-points"+yy)
                .attr("opacity", 0.1);
              grid.select(".wave-y"+yy)
                .attr("opacity", 0.1);
            }
          }

          grid.select(".line-y0")
            .transition()
            .duration(1000)
            .on("start", drawWave);
        }

        drawAnimation();

      } 

     
      drawCFLWaveAnimation();

