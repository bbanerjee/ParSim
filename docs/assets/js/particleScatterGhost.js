  var particleScatterGhost = (function() {

    var canvas = document.querySelector("#particle-scatter-ghost");
    var context = canvas.getContext("2d");
    var width = canvas.width;
    var height = canvas.height;

    var xscale = d3.scaleLinear()
                   .domain([-1, 12])
                   .range([0, width]);
    var yscale = d3.scaleLinear()
                   .domain([-1, 12])
                   .range([height, 0]);

    var line = d3.line()
                 .x(function(d) {return xscale(d.x);})
                 .y(function(d) {return yscale(d.y);});

    var initSquare = [{x: 0, y: 0},
                      {x: 3, y: 0},
                      {x: 3, y: 2},
                      {x: 0, y: 3}];

    var particles = [];
    var patches = [];
    var patchMaps = [];
    var patchColors = colorbrewer.Set1[9];

    function initialize() {
      createParticles(50);
      createPatches();
      locateParticlesInPatch();
      moveParticlesToPatches();
      createGhostRegions();
    }

    // Create particles
    function createParticles(num) {
      let count = 0;
      while (count < num) {
        let xCen = d3.randomUniform(0.1, 2.9)();
        let yCen = d3.randomUniform(0.1, 2.9)();
        if (!isOverlapping(xCen, yCen, 0.1)) {
          particles.push({x: xCen, y: yCen, moved: false});
          count++;
        }
      }
    }
        
    // Check if particle overlaps with existing
    function isOverlapping(xCen, yCen, tol) {
      for (let part of particles) {
        let distx = xCen - part.x;
        let disty = yCen - part.y;
        distSq = distx*distx + disty*disty;
        if (distSq < tol) {
          return true;
        }
      }
      return false;
    }

    // Create patches
    function createPatches() {

      let xShift = 0;
      let yShift = 0;
      let numPatch = 3;
      let patchWidth = 3; 
      let buffer = 1.2;
      let patchID = 0;
      for (let ii = 0; ii < numPatch; ii++) {
        let xmap = ii*3/3;
        let xloc = xShift + ii*(patchWidth+buffer);
        let patchXLoc = "inside";
        if (ii === 0) {
          patchXLoc = "xminus";
        } else if (ii === numPatch-1) {
          patchXLoc = "xplus";
        }
        for (let jj = 0; jj < numPatch; jj++) {
          let yloc = yShift + jj*(patchWidth+buffer);
          let ymap = jj*3/3;
          let patchYLoc = "inside";
          if (jj === 0) {
            patchYLoc = "yminus";
          } else if (jj === numPatch-1) {
            patchYLoc = "yplus";
          }
          patches.push([{x: xloc, y: yloc},
                        {x: xloc + patchWidth, y: yloc},
                        {x: xloc + patchWidth, y: yloc + patchWidth},
                        {x: xloc, y: yloc + patchWidth},
                        {x: patchXLoc, y:patchYLoc}]);
          patchMaps.push([{x: xmap, y: ymap},
                          {x: xmap + 3/3, y: ymap},
                          {x: xmap + 3/3, y: ymap + 3/3},
                          {x: xmap, y: ymap + 3/3}]);
        }
      }
    }

    // Locate particles in the patchmap
    function locateParticlesInPatch() {
      for (let part of particles) {
        let patchCount = 0;
        for (let patch of patchMaps) {
          if (inside(part, patch)) {
            let px = part.x;
            let py = part.y;
            let xmin = patch[0].x;
            let xmax = patch[1].x;
            let ymin = patch[0].y;
            let ymax = patch[3].y;
            part.patch = patchCount;
            part.tx = (px - xmin)/(xmax - xmin);
            part.ty = (py - ymin)/(ymax - ymin);
          }
          patchCount++;
        }
      }
    }

    // Find if particle is inside patch
    function inside(particle, patch) {
      if (particle.x > patch[0].x && particle.x < patch[2].x && 
          particle.y > patch[0].y && particle.y < patch[2].y) {
        return true;
      }
      return false;
    }

    // Move particles
    function moveParticlesToPatches() {

      // Compute patch locations
      for (let part of particles) {
        let patchID = 0;
        for (let patch of patches) {
          if (part.patch === patchID) {
            part.xnew = (1 - part.tx)*patch[0].x + part.tx*patch[1].x; 
            part.ynew = (1 - part.ty)*patch[0].y + part.ty*patch[3].y; 
          }
          patchID++;
        }
      }
    }

    // Create ghost regions for each patch
    function createGhostRegions() {
      let ghostWidth = 0.5;
      let patchID = 0;
      for (let patch of patches) {
        
        patch.xminus = [{x: patch[0].x, y: patch[0].y},
                        {x: patch[0].x + ghostWidth, y: patch[0].y},
                        {x: patch[0].x + ghostWidth, y: patch[3].y},
                        {x: patch[0].x, y: patch[3].y}];
        patch.xplus = [{x: patch[1].x - ghostWidth, y: patch[1].y},
                       {x: patch[1].x, y: patch[0].y},
                       {x: patch[1].x, y: patch[3].y},
                       {x: patch[1].x - ghostWidth, y: patch[3].y}];
        patch.yminus = [{x: patch[0].x, y: patch[0].y},
                        {x: patch[1].x, y: patch[0].y},
                        {x: patch[1].x, y: patch[0].y + ghostWidth},
                        {x: patch[0].x, y: patch[0].y + ghostWidth}];
        patch.yplus = [{x: patch[0].x, y: patch[3].y},
                       {x: patch[1].x, y: patch[3].y},
                       {x: patch[1].x, y: patch[3].y - ghostWidth},
                       {x: patch[0].x, y: patch[3].y - ghostWidth}];
        patchID++;
      }
    }

    // Find if particle is inside patch
    function inside(particle, patch) {
      if (particle.x > patch[0].x && particle.x < patch[2].x && 
          particle.y > patch[0].y && particle.y < patch[2].y) {
        return true;
      }
      return false;
    }

    // Move particles
    function moveParticlesToPatches() {

      // Compute patch locations
      for (let part of particles) {
        let patchID = 0;
        for (let patch of patches) {
          if (part.patch === patchID) {
            part.xnew = (1 - part.tx)*patch[0].x + part.tx*patch[1].x; 
            part.ynew = (1 - part.ty)*patch[0].y + part.ty*patch[3].y; 
          }
          patchID++;
        }
      }
    }


    function doAnimation() {
      const maxTime = 1500;
      const ease = d3.easeCubic;
      let simulation = d3.timer(
        function(elapsed) {
          const t = Math.min(1, ease(elapsed/maxTime));
   
          particles.forEach(particle => {
            let xpos = (1 - t)*particle.x + t*particle.xnew;
            let ypos = (1 - t)*particle.y + t*particle.ynew;
            particle.xpos = xpos;
            particle.ypos = ypos;
          });

          drawMovedParticles();
          drawGhostRegions();

          if (t === 1) {
            simulation.stop();
          }
        });
    }

    function reset() {
      particles = [];
      patches = [];
      patchMaps = [];
    }

    function drawBox(curSquare, dashed) {
      var coords = [];
      var index = 0;
      for (var vertex of curSquare) {
        coords[index] = {x : xscale(vertex.x),
                         y : yscale(vertex.y)};
        index++;
      }
      context.beginPath();
      if (dashed) {
        context.setLineDash([5]);
      }
      context.moveTo(coords[0].x, coords[0].y);
      context.lineTo(coords[1].x, coords[1].y);
      context.lineTo(coords[2].x, coords[2].y);
      context.lineTo(coords[3].x, coords[3].y);
      context.lineTo(coords[0].x, coords[0].y);
      context.stroke();
      context.fill();
      if (dashed) {
        context.setLineDash([0]);
      }
    }

    // Draw patches
    function drawPatches() {
      let patchID = 0;
      for (let patch of patches) {
        context.save();
        context.lineWidth = 2;
        //if (patchID === 0) {
        //  context.fillStyle = "#87ceeb";
        //} else {
        //  context.fillStyle = "#9acd32";
        //}
        context.fillStyle = patchColors[patchID];
        context.globalAlpha = 0.3;
        drawBox(patch, false);

        let label = "P" + patchID;
        context.font = "16px Arial";
        context.fillStyle = "#aa0902";
        context.fillText(label, xscale(patch[0].x)+5, yscale(patch[0].y)-5);
        context.stroke();
        context.restore();
        patchID++;
      }
    }

    // Draw moved particles
    function drawMovedParticles() {

      context.clearRect(0, 0, width, height);

      drawPatches();

      for (let part of particles) {
        context.beginPath()
        context.arc(xscale(part.xpos), yscale(part.ypos), 5, 0, 2.0*Math.PI); 
        //context.fillStyle = "#aa0902";
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();
      }
    }

    function drawGhostRegions() {
      context.save();
      context.lineWidth = 2;
      context.globalAlpha = 0.4;
      let patchID = 0;
      for (let patch of patches) {
        context.fillStyle = patchColors[patchID];
        if (patch[4].x != "xminus") {
          drawBox(patch.xminus, true);
        }
        if (patch[4].x != "xplus") {
          drawBox(patch.xplus, true);
        }
        if (patch[4].y != "yminus") {
          drawBox(patch.yminus, true);
        }
        if (patch[4].y != "yplus") {
          drawBox(patch.yplus, true);
        }
        patchID++;
      }
      context.restore();
    }

    return {
      run : function() {
        initialize();
        doAnimation();
      },
      restartAnimation : function() {
        reset();
        initialize();
        doAnimation();
      }
    };
    
  })();

  particleScatterGhost.run();
