      "use strict";
      function computeYieldSurface(iterData, numPts) {

        // Get the moduli
        let K = iterData['K'];
        let G = iterData['G'];

        // Get yield parameters
        let PEAKI1 = iterData['PEAKI1'];
        let FSLOPE = iterData['FSLOPE'];
        let STREN  = iterData['STREN'];
        let YSLOPE = iterData['YSLOPE'];
        let CR     = iterData['CR'];

        // Set up constants
        let a1 = STREN;
        let a2 = (FSLOPE-YSLOPE)/(STREN-YSLOPE*PEAKI1);
        let a3 = (STREN-YSLOPE*PEAKI1)*Math.exp(-a2*PEAKI1);
        let a4 = YSLOPE;

        // Compute kappa
        let X_eff = iterData['X'] + 3.0*iterData['pbar_w'];
        let kappa = PEAKI1 - CR*(PEAKI1 - X_eff);

        // Create an array of I1_eff values
        let I1eff_min = X_eff;
        let I1eff_max = PEAKI1;

        let rad = 0.5*(PEAKI1 - X_eff);
        let cen = 0.5*(PEAKI1 + X_eff);
        let theta_max = Math.acos(Math.max((I1eff_min - cen)/rad, -1.0));
        let theta_min = Math.acos(Math.min((I1eff_max - cen)/rad, 1.0));
        let theta_inc = (theta_max-theta_min)/numPts;
        let theta_vec = d3.range(theta_min, theta_max+theta_inc, theta_inc);
        let I1_list = theta_vec.map(function(theta) {return cen + rad*Math.cos(theta);});
        let I1_eff_list = I1_list.map(function(I1) {return Math.max(I1, X_eff);});

        // Create an array of sqrtJ2 values
        let sqrtJ2_list = I1_eff_list.map(
          function(I1_eff) {
            // Compute F_f
            let Ff = a1 - a3*Math.exp(a2*I1_eff) - a4*(I1_eff);

            // Compute Fc
            let Fc_sq = 1.0;
            if ((I1_eff < kappa) && (X_eff <= I1_eff)) {
              let ratio = (kappa - I1_eff)/(kappa - X_eff);
              Fc_sq = 1.0 - ratio*ratio;
            }

            // Compute sqrt(J2)
            let sqrtJ2 = Ff*Math.sqrt(Fc_sq);
            return(sqrtJ2);
          });

        // Compute z and r'
        let z_list = I1_eff_list.map(function(I1_eff) {return I1_eff/Math.sqrt(3.0);});
        let rprime_list = sqrtJ2_list.map(function(sqrtJ2) {return sqrtJ2*Math.sqrt(3.0*K/G);});

        // Create data frames
        let zrprime_data = d3.zip(z_list, rprime_list);
        return zrprime_data;
      }

      function drawYieldSurface(data) {

        // Create svg
        var svgsize = {x: 600, y: 600};
        var margin = {top: 50, right: 50, bottom: 50, left: 50};
        var width = svgsize.x - margin.left - margin.right;
        var height = svgsize.y - margin.top - margin.bottom;

        var svg = d3.select(".yield-surf-canvas")
                    .append("svg")
                      .attr("width", svgsize.x)
                      .attr("height", svgsize.y);
        var chart = svg.append("g")
                       .attr("class", "chart")
                       .attr("transform", 
                             "translate(" + margin.left + "," + margin.top + ")");

        // Process the data
        let numPts = data['num_points'];
        let numIter = data.data.length;

        // Compute yield surface in z and r' from first element of data
        let iterData = data.data[0];
        let zrprime = computeYieldSurface(iterData, numPts);

        // Get the coordinates of the trial stress
        let ztrial = iterData.z_r_pt[0];
        let rprimetrial = iterData.z_r_pt[1];

        // Compute min max of domain
        let zmin = d3.min(zrprime, function(d) {return d[0];})
        let zmax = d3.max(zrprime, function(d) {return d[0];})
        let rprimemin = d3.min(zrprime, function(d) {return d[1];})
        let rprimemax = d3.max(zrprime, function(d) {return d[1];})
        zmin = Math.min(zmin, ztrial);
        zmax = Math.max(zmax, ztrial);
        rprimemin = Math.min(rprimemin, rprimetrial);
        rprimemax = Math.max(rprimemax, rprimetrial);

        // Create scaling functions
        var xscale = d3.scaleLinear().domain([zmin, zmax]).rangeRound([0, width]);
        var yscale = d3.scaleLinear().domain([rprimemin, rprimemax]).rangeRound([height, 0]);
        var xscaleAxis = d3.scaleLinear().domain([zmin*1.0e-6, zmax*1.0e-6]).rangeRound([0, width]);
        var yscaleAxis = d3.scaleLinear().domain([rprimemin*1.0e-6, rprimemax*1.0e-6]).rangeRound([height, 0]);

        // Create polyline generation function
        var line = d3.line()
                     .x(function(d) { return xscale(d[0]); })
                     .y(function(d) { return yscale(d[1]); });

        // Create bottom axis
        chart.append("g")
               .attr("transform", "translate(0," + height + ")")
             .call(d3.axisBottom(xscaleAxis))
             .append("text")
                .attr("fill", "#000")
               .attr("transform",
                     "translate(" + (width/2) + " ," + 
                             (margin.top - 20) + ")")
               .attr("text-anchor", "end")
               .text("z (MPa)"); 

        // Create left axis
        chart.append("g")
             .call(d3.axisLeft(yscaleAxis))
           .append("text")
             .attr("fill", "#000")
             .attr("transform", "rotate(-90)")
             .attr("y", 6)
             .attr("dy", "0.71em")
             .attr("text-anchor", "end")
             .text("r' (MPa)"); 

        // Plot yield surface
        chart.append("path")
             .attr("class", "yield_surf")
             .datum(zrprime)
             .attr("fill", "none")    
             .attr("stroke", "steelblue")
             .attr("stroke-linejoin", "round")
             .attr("stroke-linecap", "round")
             .attr("stroke-width", 1.5)
             .attr("d", line);

        // Loop through the iterations
        for (let iter = 0; iter < numIter; iter++) {

          iterData = data.data[iter];

          // Create group for each iteration
          let iterGroup = chart.append("g")
                               .attr("class", "iteration"+iter)
                               .attr("opacity", 0);

          // Create group for the yield surface points + circles
          let yieldSurfGroup = iterGroup.append("g") 
                                        .attr("class", "approx_yield_surf"+iter);

          // Get the cordinates of polyline approximation
          let zpoly = iterData.z_r_yield_z;
          let rprimepoly = iterData.z_r_yield_r;
          let zrprimepoly = d3.zip(zpoly, rprimepoly);

          yieldSurfGroup.append("path")
             .datum(zrprimepoly)
             .attr("fill", "none")    
             .attr("stroke", "red")
             .attr("stroke-linejoin", "round")
             .attr("stroke-linecap", "round")
             .attr("stroke-width", 1.5)
             .attr("d", line);

          yieldSurfGroup.selectAll(".approx_yield_surf_circle"+iter)
            .data(zrprimepoly)
            .enter()
            .append("circle")
            .attr("r", 2)
            .attr("cx", function(d) {return xscale(d[0]);})
            .attr("cy", function(d) {return yscale(d[1]);})
            .attr("fill", "blue")
            .attr("stroke", "black");

          // Create group for the projection line and points
          let projLineGroup = iterGroup.append("g") 
                                       .attr("class", "projection_line"+iter);

          // Get the coordinates of the trial stress
          let ztrial = iterData.z_r_pt[0];
          let rprimetrial = iterData.z_r_pt[1];

          // Get the coordinates of the closest point stress
          let zclosest = iterData.z_r_closest[0];
          let rprimeclosest = iterData.z_r_closest[1];

          // Set up projection line
          let zrprimeproj = [[ztrial, rprimetrial],[zclosest, rprimeclosest]];

          projLineGroup.append("path")
            .datum(zrprimeproj)
            .attr("fill", "none")    
            .attr("stroke", "green")
            .attr("stroke-linejoin", "round")
            .attr("stroke-linecap", "round")
            .attr("stroke-width", 1.5)
            .attr("d", line);

          projLineGroup.selectAll(".projection_line_point"+iter)
            .data(zrprimeproj)
            .enter()
            .append("path")
            .attr("class", "projection_line_point"+iter)
            .attr("d", d3.symbol().type(d3.symbolTriangle))
            .attr("transform", 
              function(d) {return "translate(" + xscale(d[0]) + "," + yscale(d[1]) + ")";});

        }

        function drawAnimation() {
          var iterationID = 0;

          function drawNextSurf() {
            d3.active(this)
                .attr("opacity", 1)
              .transition()
                .duration(1000)
                .attr("opacity", 0);
            iterationID++;
            if (iterationID === numIter) {
              iterationID = 0;
            }
            chart.select(".iteration"+iterationID)
                 .transition()
                   .delay(1000)
                   .duration(1000)
                 .on("start", drawNextSurf);            
          }

          chart.select(".iteration0")
               .transition()
                 .duration(1000)
               .on("start", drawNextSurf);
        }

        drawAnimation();

      } 
     
