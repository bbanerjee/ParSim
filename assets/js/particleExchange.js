  var particleExchange = (function () {

    var canvas = document.querySelector("#particle-exchange-ghost");
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

    // Patch display variables
    var xShift = 0;
    var yShift = 0;
    var numPatch = 3;
    var patchWidth = 3; 
    var buffer = 1.2;
    var ghostWidth = 0.5;

    function initialize() {
      createParticles(50);
      createPatches();
      locateParticlesInPatch();
      moveParticlesToPatches();
      createGhostRegions();
      addParticlesToGhostRegions();
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
    function inside(particle, patch, moved = false) {
      if (!moved) {
        if (particle.x > patch[0].x && particle.x < patch[2].x && 
            particle.y > patch[0].y && particle.y < patch[2].y) {
          return true;
        }
      } else {
        if (particle.xnew > patch[0].x && particle.xnew < patch[2].x && 
            particle.ynew > patch[0].y && particle.ynew < patch[2].y) {
          return true;
        }
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
        patch.yplus = [{x: patch[0].x, y: patch[3].y - ghostWidth},
                       {x: patch[1].x, y: patch[3].y - ghostWidth},
                       {x: patch[1].x, y: patch[3].y},
                       {x: patch[0].x, y: patch[3].y}];
        patch.xminusyminus = [{x: patch[0].x, y: patch[0].y},
                              {x: patch[0].x + ghostWidth, y: patch[0].y},
                              {x: patch[0].x + ghostWidth, y: patch[0].y + ghostWidth},
                              {x: patch[0].x, y: patch[0].y + ghostWidth}];
        patch.xminusyplus = [{x: patch[0].x, y: patch[3].y - ghostWidth},
                             {x: patch[0].x + ghostWidth, y: patch[3].y - ghostWidth},
                             {x: patch[0].x + ghostWidth, y: patch[3].y},
                             {x: patch[0].x, y: patch[3].y}];
        patch.xplusyminus = [{x: patch[1].x - ghostWidth, y: patch[0].y},
                             {x: patch[1].x, y: patch[0].y},
                             {x: patch[1].x, y: patch[0].y + ghostWidth},
                             {x: patch[1].x - ghostWidth, y: patch[0].y + ghostWidth}];
        patch.xplusyplus = [{x: patch[1].x - ghostWidth, y: patch[3].y - ghostWidth},
                            {x: patch[1].x, y: patch[3].y - ghostWidth},
                            {x: patch[1].x, y: patch[3].y},
                            {x: patch[1].x - ghostWidth, y: patch[3].y}];

        patchID++;
      }
    }

    function addParticlesToGhostRegions() {
      let moved = true;
      for (let patch of patches) {
        patch.xminus.particles = [];
        patch.xplus.particles = [];
        patch.yminus.particles = [];
        patch.yplus.particles = [];
        patch.xminusyminus.particles = [];
        patch.xminusyplus.particles = [];
        patch.xplusyminus.particles = [];
        patch.xplusyplus.particles = [];
        for (let particle of particles) {
          if (inside(particle, patch.xminus, moved)) {
            patch.xminus.particles.push(particle);
          }
          if (inside(particle, patch.xplus, moved)) {
            patch.xplus.particles.push(particle);
          }
          if (inside(particle, patch.yminus, moved)) {
            patch.yminus.particles.push(particle);
          }
          if (inside(particle, patch.yplus, moved)) {
            patch.yplus.particles.push(particle);
          }
          if (inside(particle, patch.xminusyminus, moved)) {
            patch.xminusyminus.particles.push(particle);
          }
          if (inside(particle, patch.xminusyplus, moved)) {
            patch.xminusyplus.particles.push(particle);
          }
          if (inside(particle, patch.xplusyminus, moved)) {
            patch.xplusyminus.particles.push(particle);
          }
          if (inside(particle, patch.xplusyplus, moved)) {
            patch.xplusyplus.particles.push(particle);
          }
        }
      }
    }

    function moveGhostRegions(tt) {
      for (let patch of patches) {
        if (patch[4].x != "xminus") {
          moveXMinus(patch, tt);
        }
        if (patch[4].x != "xplus") {
          moveXPlus(patch, tt);
        }
        if (patch[4].y != "yminus") {
          moveYMinus(patch, tt);
        }
        if (patch[4].y != "yplus") {
          moveYPlus(patch, tt);
        }
        if (patch[4].x != "xminus" && patch[4].y != "yminus") {
          moveXMinusYMinus(patch, tt);
        }
        if (patch[4].x != "xminus" && patch[4].y != "yplus") {
          moveXMinusYPlus(patch, tt);
        }
        if (patch[4].x != "xplus" && patch[4].y != "yminus") {
          moveXPlusYMinus(patch, tt);
        }
        if (patch[4].x != "xplus" && patch[4].y != "yplus") {
          moveXPlusYPlus(patch, tt);
        }
      }
    }

    function moveXMinus(patch, tt) {
      for (let ii = 0; ii < 4; ii++) {
        patch.xminus[ii].xghost = patch.xminus[ii].x - tt*buffer;
        patch.xminus[ii].yghost = patch.xminus[ii].y;
      }
      for (let particle of patch.xminus.particles) {
        particle.xghostxm = particle.xnew - tt*buffer;
        particle.yghostxm = particle.ynew;
      }
    }

    function moveXPlus(patch, tt) {
      for (let ii = 0; ii < 4; ii++) {
        patch.xplus[ii].xghost = patch.xplus[ii].x + tt*buffer;
        patch.xplus[ii].yghost = patch.xplus[ii].y;
      }
      for (let particle of patch.xplus.particles) {
        particle.xghostxp = particle.xnew + tt*buffer;
        particle.yghostxp = particle.ynew;
      }
    }

    function moveYMinus(patch, tt) {
      for (let ii = 0; ii < 4; ii++) {
        patch.yminus[ii].xghost = patch.yminus[ii].x;
        patch.yminus[ii].yghost = patch.yminus[ii].y - tt*buffer;
      }
      for (let particle of patch.yminus.particles) {
        particle.xghostym = particle.xnew;
        particle.yghostym = particle.ynew - tt*buffer;
      }
    }

    function moveYPlus(patch, tt) {
      for (let ii = 0; ii < 4; ii++) {
        patch.yplus[ii].xghost = patch.yplus[ii].x;
        patch.yplus[ii].yghost = patch.yplus[ii].y + tt*buffer;
      }
      for (let particle of patch.yplus.particles) {
        particle.xghostyp = particle.xnew;
        particle.yghostyp = particle.ynew + tt*buffer;
      }
    }

    function moveXMinusYMinus(patch, tt) {
      for (let ii = 0; ii < 4; ii++) {
        patch.xminusyminus[ii].xghost = patch.xminusyminus[ii].x - tt*buffer;
        patch.xminusyminus[ii].yghost = patch.xminusyminus[ii].y - tt*buffer;
      }
      for (let particle of patch.xminusyminus.particles) {
        particle.xghostxmym = particle.xnew - tt*buffer;
        particle.yghostxmym = particle.ynew - tt*buffer;
      }
    }

    function moveXMinusYPlus(patch, tt) {
      for (let ii = 0; ii < 4; ii++) {
        patch.xminusyplus[ii].xghost = patch.xminusyplus[ii].x - tt*buffer;
        patch.xminusyplus[ii].yghost = patch.xminusyplus[ii].y + tt*buffer;
      }
      for (let particle of patch.xminusyplus.particles) {
        particle.xghostxmyp = particle.xnew - tt*buffer;
        particle.yghostxmyp = particle.ynew + tt*buffer;
      }
    }

    function moveXPlusYMinus(patch, tt) {
      for (let ii = 0; ii < 4; ii++) {
        patch.xplusyminus[ii].xghost = patch.xplusyminus[ii].x + tt*buffer;
        patch.xplusyminus[ii].yghost = patch.xplusyminus[ii].y - tt*buffer;
      }
      for (let particle of patch.xplusyminus.particles) {
        particle.xghostxpym = particle.xnew + tt*buffer;
        particle.yghostxpym = particle.ynew - tt*buffer;
      }
    }

    function moveXPlusYPlus(patch, tt) {
      for (let ii = 0; ii < 4; ii++) {
        patch.xplusyplus[ii].xghost = patch.xplusyplus[ii].x + tt*buffer;
        patch.xplusyplus[ii].yghost = patch.xplusyplus[ii].y + tt*buffer;
      }
      for (let particle of patch.xplusyplus.particles) {
        particle.xghostxpyp = particle.xnew + tt*buffer;
        particle.yghostxpyp = particle.ynew + tt*buffer;
      }
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
      const maxTime = 3000;
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

          moveGhostRegions(t);
          drawMovedParticles();
          drawMovedGhostRegions();
          drawMovedGhostParticles();
          drawMovedGhostCornerParticles();

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

    function drawBox(curSquare, dashed, type="none") {
      var coords = [];
      var index = 0;
      if (type === "ghost") {
        for (var vertex of curSquare) {
          coords[index] = {x : xscale(vertex.xghost),
                           y : yscale(vertex.yghost)};
          index++;
        }
      } else {
        for (var vertex of curSquare) {
          coords[index] = {x : xscale(vertex.x),
                           y : yscale(vertex.y)};
          index++;
        }
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
        context.fillStyle = patchColors[patchID];
        //if (patchID === 0) {
        //  context.fillStyle = "#87ceeb";
        //} else {
        //  context.fillStyle = "#9acd32";
        //}
        context.globalAlpha = 0.2;
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
        context.arc(xscale(part.xnew), yscale(part.ynew), 5, 0, 2.0*Math.PI); 
        //context.fillStyle = "#aa0902";
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();
      }
    }

    function drawMovedGhostRegions() {
      context.save();
      context.lineWidth = 2;
      context.globalAlpha = 0.4;
      let patchID = 0;
      for (let patch of patches) {
        context.fillStyle = patchColors[patchID];
        if (patch[4].x != "xminus") {
          drawBox(patch.xminus, true, "ghost");
        }
        if (patch[4].x != "xplus") {
          drawBox(patch.xplus, true, "ghost");
        }
        if (patch[4].y != "yminus") {
          drawBox(patch.yminus, true, "ghost");
        }
        if (patch[4].y != "yplus") {
          drawBox(patch.yplus, true, "ghost");
        }
        if (patch[4].x != "xminus" && patch[4].x != "yminus") {
          drawBox(patch.xminusyminus, true, "ghost");
        }
        if (patch[4].x != "xminus" && patch[4].x != "yplus") {
          drawBox(patch.xminusyplus, true, "ghost");
        }
        if (patch[4].x != "xplus" && patch[4].x != "yminus") {
          drawBox(patch.xplusyminus, true, "ghost");
        }
        if (patch[4].x != "xplus" && patch[4].x != "yplus") {
          drawBox(patch.xplusyplus, true, "ghost");
        }
        patchID++;
      }
      context.restore();
    }

    function drawMovedGhostParticles() {

      for (let part of particles) {
        context.beginPath()
        context.arc(xscale(part.xghostxm), yscale(part.yghostxm), 5, 0, 2.0*Math.PI); 
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();

        context.beginPath()
        context.arc(xscale(part.xghostxp), yscale(part.yghostxp), 5, 0, 2.0*Math.PI); 
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();

        context.beginPath()
        context.arc(xscale(part.xghostym), yscale(part.yghostym), 5, 0, 2.0*Math.PI); 
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();

        context.beginPath()
        context.arc(xscale(part.xghostyp), yscale(part.yghostyp), 5, 0, 2.0*Math.PI); 
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();
      }
    }

    function drawMovedGhostCornerParticles() {

      for (let part of particles) {
        context.beginPath()
        context.arc(xscale(part.xghostxmym), yscale(part.yghostxmym), 5, 0, 2.0*Math.PI); 
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();

        context.beginPath()
        context.arc(xscale(part.xghostxmyp), yscale(part.yghostxmyp), 5, 0, 2.0*Math.PI); 
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();

        context.beginPath()
        context.arc(xscale(part.xghostxpym), yscale(part.yghostxpym), 5, 0, 2.0*Math.PI); 
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();

        context.beginPath()
        context.arc(xscale(part.xghostxpyp), yscale(part.yghostxpyp), 5, 0, 2.0*Math.PI); 
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();
      }
    }

    return {
      run : function () {
        initialize();
        doAnimation();
      },
      restartAnimation: function () {
        reset();
        initialize();
        doAnimation();
      }
    };
  })();

  particleExchange.run();
    


