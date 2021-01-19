
      function drawCFLDomainAnimation() {

        // Create svg
        let svgsize = {x: 500, y: 300};
        let margin = {top: 50, right: 50, bottom: 50, left: 50};
        let width = svgsize.x - margin.left - margin.right;
        let height = svgsize.y - margin.top - margin.bottom;

        let svg = d3.select(".cfl-domain-animation")
                    .append("svg")
                    .attr("width", svgsize.x)
                    .attr("height", svgsize.y);
        let grid = svg.append("g")
                     .attr("class", "grid")
                     .attr("transform", 
                           "translate(" + margin.left + "," + margin.top + ")");
        let slider = svg.append("g")
                       .attr("class", "slider")
                       .attr("transform", 
                             "translate(" + margin.left + "," + (margin.top/2 - 10) + ")");

        // Create scaling functions
        let xmin = -10.0;
        let xmax =  10.0;
        let ymin =  0.0;
        let ymax =  5.0;
        let xscale = d3.scaleLinear().domain([xmin, xmax]).rangeRound([0, width]);
        let yscale = d3.scaleLinear().domain([ymin, ymax]).rangeRound([height, 0]);
        let kmin = 0.5;
        let kmax = 2.0;
        let kscale = d3.scaleLinear().domain([kmin, kmax]).rangeRound([0, width]).clamp(true);

        // Create polyline generation function
        let line = d3.line()
                     .x(function(d) { return xscale(d[0]); })
                     .y(function(d) { return yscale(d[1]); });

        // Create slider track
        let kval = 1.0;
        let kchanged = false;
        let kticks = kscale.ticks(4);
        let kindex = 1;
        slider.append("line")
          .attr("class", "track")
          .attr("x1", kscale.range()[0])
          .attr("x2", kscale.range()[1])
          .select(function() {
                    return this.parentNode.appendChild(this.cloneNode(true));
                  })
          .attr("class", "track-inset")
          .select(function() {
                    return this.parentNode.appendChild(this.cloneNode(true));
                  })
          .attr("class", "track-overlay")
          .call(d3.drag()
                  .on("start.interrupt", 
                       function() {resetAll(); slider.interrupt(); })
                  .on("start drag", 
                      function() {
                        kval = kscale.invert(d3.event.x);
                        let tt = kticks.map(function(k) {return (kval-k)/0.5;});
                        let ttlocLo = tt.map(function(t) {return (t >= 0.0 && t < 0.5);});
                        let ttlocHi = tt.map(function(t) {return (t >= 0.5 && t <= 1.0);});
                        let ttlocLoTrue = ttlocLo.reduce(function(a, b) {return (a | b);});
                        let ttlocHiTrue = ttlocHi.reduce(function(a, b) {return (a | b);});
                        if (ttlocLoTrue) {
                          for (let ii = 0; ii < 4; ii++) {
                            if (ttlocLo[ii] == true) {
                              kindex = ii;
                              break;
                            }
                          }
                        } else if (ttlocHiTrue) {
                          for (let ii = 0; ii < 4; ii++) {
                            if (ttlocHi[ii] == true) {
                              kindex = ii+1;
                              break;
                            }
                          }
                        }
                        kval = kticks[kindex];
                        handle.attr("cx", kscale(kval));
                      }));

        // Add the slider labels
        slider.insert("g", ".track-overlay")
          .attr("class", "ticks")
          .attr("transform", "translate(0," + 18 + ")")
          .selectAll("text")
          .data(kticks)
          .enter().append("text")
          .attr("x", kscale)
          .attr("text-anchor", "middle")
          .text(function(d) { return d3.format(".2n")(d); });

        slider.insert("g", ".track-text")
          .attr("class", "slider-name")
          .attr("transform", "translate(-20," + "5" + ")")
          .append("text")
          .attr("x", kscale(kticks[0]))
          .attr("text-anchor", "middle")
          .text("k");
                       
        // Add the slider handle
        let handle = slider.insert("circle", ".track-overlay")
                       .attr("class", "handle")
                       .attr("r", 9)
                       .attr("cx", kscale(1.0));
       
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

        
        // Get the x and y increments for each case
        let yinc = 1.0;
        let xinc = kticks.map(function(k) {return yinc*k;});

        // Plot horizontal lines
        for (let yCoord = ymin; yCoord < ymax+1; yCoord++) {
          grid.append("g")
            .append("path")
            .attr("class", "line-y"+yCoord)
            .datum([[xmin, yCoord],[xmax, yCoord]])
            .attr("fill", "none")    
            .attr("stroke", "gray")
            .attr("stroke-width", 0.5)
            .attr("opacity", 0)
            .attr("d", line);
        }

        // Create geometries for each of the values of k
        for (let kk = 0; kk < 4; kk++) {

          // Plot grid points
          let gridPoints = [];
          for (let yCoord = ymin; yCoord < ymax+1; yCoord++) {
            let gridData = [];
            let index = 0;
            for (let xCoord = 0; xCoord > xmin; xCoord-=xinc[kk]) {
              gridData[index] = [xCoord, yCoord];
              index++;
            }
            for (let xCoord = 0; xCoord < xmax; xCoord+=xinc[kk]) {
              gridData[index] = [xCoord, yCoord];
              index++;
            }
            gridPoints[yCoord] = grid.append("g")
              .attr("class", "grid-points"+yCoord+"-"+kk)
              .attr("opacity", 0)
              .selectAll("grid-point-y"+yCoord+"-"+kk)
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
            waveFrontData = [[(- yCoord*kticks[kk]), ymin],
                             [( yCoord*kticks[kk]), ymin],
                             [0, yCoord]]; 
            waveFront[yCoord] = grid.append("g")
              .append("path")
              .attr("class", "wave-y"+yCoord+"-"+kk)
              .datum(waveFrontData)
              .attr("fill", "pink")    
              .attr("stroke", "gray")
              .attr("stroke-width", 0.5)
              .attr("opacity", 0)
              .attr("d", line);
          }

          // Plot wave front mathematical
          let waveFrontMath = [];
          for (let yCoord = ymin; yCoord < ymax+1; yCoord++) {
            let waveFrontData = [];
            waveFrontData = [[(- yCoord), ymin],
                             [( yCoord), ymin],
                             [0, yCoord],
                             [(- yCoord), ymin]]; 
            waveFrontMath[yCoord] = grid.append("g")
              .append("path")
              .attr("class", "wave-math-y"+yCoord+"-"+kk)
              .datum(waveFrontData)
              .attr("fill", "none")    
              .attr("stroke", "green")
              .attr("stroke-width", 1.0)
              .attr("opacity", 1)
              .attr("d", line);
          }
                
        }

        let currY = ymin;

        function resetAll() {
          currY = ymin;
          for (let kk = 0; kk < 4; kk++) {
            for (let yy = ymin; yy < ymax+1; yy++) {
              grid.select(".line-y"+yy)
                .attr("opacity", 0.1);
              grid.select(".grid-points"+yy+"-"+kk)
                .attr("opacity", 0.0);
              grid.select(".wave-y"+yy+"-"+kk)
                .attr("opacity", 0.0);
            }
          }
        }

        function drawAnimation() {

          currY = ymin;
          function drawWave() {
            d3.active(this)
              .attr("opacity", 1);
            grid.selectAll(".grid-points"+currY+"-"+kindex)
              .attr("opacity", 0.8);
            grid.select(".wave-y"+currY+"-"+kindex)
              .attr("opacity", 0.6);
            currY++;
            if (currY > ymax) {
              currY = ymin;
              resetWave();
            }
            grid.select(".line-y"+currY)
              .transition()
              .delay(1000)
              .duration(1000)
              .on("start", drawWave);
          }

          function resetWave() {
            for (let yy = ymin; yy < ymax+1; yy++) {
              grid.select(".line-y"+yy)
                .attr("opacity", 0.1);
              grid.selectAll(".grid-points"+yy+"-"+kindex)
                .attr("opacity", 0.0);
              grid.select(".wave-y"+yy+"-"+kindex)
                .attr("opacity", 0.0);
            }
          }

          grid.select(".line-y0")
            .transition()
            .duration(1000)
            .on("start", drawWave);
        }

        drawAnimation();

      } 

     
      drawCFLDomainAnimation();

