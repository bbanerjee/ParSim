  var particleMotion = (function () {

    var canvas = document.querySelector("#particle-motion");

    var context = canvas.getContext("2d");
    var width = canvas.width;
    var height = canvas.height;

    /*
    var hiddenCanvas = d3.select("#container")
      .append('canvas')
      .classed('hiddenCanvas', true)
      .attr('width', width)
      .attr('height', height);
    */

    var xscale;
    var yscale;

    var line = d3.line()
                 .x(function(d) {return xscale(d.x);})
                 .y(function(d) {return yscale(d.y);});

    var initSquare = [{x: 0, y: 0},
                      {x: 3, y: 0},
                      {x: 3, y: 2},
                      {x: 0, y: 3}];

    var particles = [];
    var patches = [];
    var patchesExtended = [];
    var patchMaps = [];
    var patchColors = colorbrewer.Set1[9];

    // Patch display variables
    var xShift = 0;
    var yShift = 0;
    var numPatch = 3;
    var patchWidth = 3; 
    var buffer;
    var ghostWidth;

    function initialize() {
      initScales();
      initGhostWidth();
      createParticles(50);
      createPatches();
      locateParticlesInPatch();
      moveParticlesToPatches();
      createGhostRegions();
      addParticlesToGhostRegionsX();
      addParticlesToGhostRegionsY();

      // Set up the expanded patches first
      context.clearRect(0, 0, width, height);
      drawPatches();
      moveGhostRegionsX(1.0);
      drawMovedParticles();
      drawMovedGhostRegionsX();
      drawMovedGhostParticlesX();
      appendParticlesToGhostRegionsY();
      moveGhostRegionsY(1.0);
      drawMovedParticles();
      drawMovedGhostRegionsX();
      drawMovedGhostParticlesX();
      drawCopiedGhostParticlesX(1.0);
      drawMovedGhostRegionsY();
      drawMovedGhostParticlesY();

      // Create the updated particle list for each patch
      createUpdatedParticleLists();
    }

    function initScales() {
      xscale = d3.scaleLinear()
                   .domain([-1, 14])
                   .range([0, width]);
      yscale = d3.scaleLinear()
                   .domain([-1, 14])
                   .range([height, 0]);
    }

    function initGhostWidth() {
      buffer = 2.0;
      ghostWidth = 0.5;
    }

    // Create particles
    function createParticles(num) {
      let seed = Math.seedrandom("12345");
      let count = 0;
      while (count < num) {
        let val = Math.random();
        let xCen = (1 - val)*0.1 + val*2.9;
        val = Math.random();
        let yCen = (1 - val)*0.1 + val*2.9;
        val = Math.random();
        let xVel = (1 - val)*(-1.0) + val*1.0;
        val = Math.random();
        let yVel =  (1 - val)*(-1.0) + val*1.0;
        if (!isOverlapping(xCen, yCen, 0.1)) {
          particles.push({x: xCen, y: yCen, moved: false, id: count, xvel: xVel, yvel: yVel});
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
          patchesExtended.push([{x: xloc-ghostWidth, y: yloc-ghostWidth},
                                {x: xloc+ghostWidth + patchWidth, y: yloc-ghostWidth},
                                {x: xloc+ghostWidth + patchWidth, y: yloc+ghostWidth + patchWidth},
                                {x: xloc-ghostWidth, y: yloc+ghostWidth + patchWidth}]);
          patchMaps.push([{x: xmap, y: ymap},
                          {x: xmap + 3/3, y: ymap},
                          {x: xmap + 3/3, y: ymap + 3/3},
                          {x: xmap, y: ymap + 3/3}]);
        }
      }
      let patchID = 0;
      for (let patch of patches) {
        patch.id = patchID;
        patch.particles = [];
        patch.xminusparticles = [];
        patch.xplusparticles = [];
        patch.yminusparticles = [];
        patch.yplusparticles = [];
        patch.yminusxminusparticles = [];
        patch.yminusxplusparticles = [];
        patch.yplusxminusparticles = [];
        patch.yplusxplusparticles = [];
        patch.totalparticles = [];
        patch.updated = [];
        patchID++;
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
            patch.particles.push(part);
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
        patch.yminus = [{x: patch[0].x - ghostWidth, y: patch[0].y},
                        {x: patch[1].x + ghostWidth, y: patch[0].y},
                        {x: patch[1].x + ghostWidth, y: patch[0].y + ghostWidth},
                        {x: patch[0].x - ghostWidth, y: patch[0].y + ghostWidth}];
        patch.yplus = [{x: patch[0].x - ghostWidth, y: patch[3].y - ghostWidth},
                       {x: patch[1].x + ghostWidth, y: patch[3].y - ghostWidth},
                       {x: patch[1].x + ghostWidth, y: patch[3].y},
                       {x: patch[0].x - ghostWidth, y: patch[3].y}];

        patchID++;
      }
    }

    function addParticlesToGhostRegionsX() {
      for (let patch of patches) {
        patch.xminus.particles = [];
        patch.xplus.particles = [];
        for (let particle of particles) {
          let pp = {x: particle.xnew, y:particle.ynew};
          if (inside(pp, patch.xminus)) {
            patch.xminusparticles.push(particle);
          }
          if (inside(pp, patch.xplus)) {
            patch.xplusparticles.push(particle);
          }
        }
      }
    }

    function addParticlesToGhostRegionsY() {
      for (let patch of patches) {
        patch.yminus.particles = [];
        patch.yplus.particles = [];
        for (let particle of particles) {
          let pp = {x: particle.xnew, y: particle.ynew};
          if (inside(pp, patch.yminus)) {
            patch.yminusparticles.push(particle);
          }
          if (inside(pp, patch.yplus)) {
            patch.yplusparticles.push(particle);
          }
        }
      }
    }

    function appendParticlesToGhostRegionsY() {
      for (let patch of patches) {
        for (let particle of particles) {
          let pp = {x: particle.xghostxm, y:particle.yghostxm};
          if (inside(pp, patch.yminus)) {
            patch.yminusxminusparticles.push(particle);
          }
          if (inside(pp, patch.yplus)) {
            patch.yplusxminusparticles.push(particle);
          }
          pp = {x: particle.xghostxp, y:particle.yghostxp};
          if (inside(pp, patch.yminus)) {
            patch.yminusxplusparticles.push(particle);
          }
          if (inside(pp, patch.yplus)) {
            patch.yplusxplusparticles.push(particle);
          }
        }
      }
    }

    function moveGhostRegionsX(tt) {
      for (let patch of patches) {
        if (patch[4].x != "xminus") {
          moveXMinus(patch, tt);
        }
        if (patch[4].x != "xplus") {
          moveXPlus(patch, tt);
        }
      }
    }

    function moveGhostRegionsY(tt) {
      for (let patch of patches) {
        if (patch[4].y != "yminus") {
          moveYMinus(patch, tt);
        }
        if (patch[4].y != "yplus") {
          moveYPlus(patch, tt);
        }
      }
    }

    function moveXMinus(patch, tt) {
      for (let ii = 0; ii < 4; ii++) {
        patch.xminus[ii].xghost = patch.xminus[ii].x - tt*buffer;
        patch.xminus[ii].yghost = patch.xminus[ii].y;
      }
      for (let particle of patch.xminusparticles) {
        particle.xghostxm = particle.xnew - tt*buffer;
        particle.yghostxm = particle.ynew;
      }
    }

    function moveXPlus(patch, tt) {
      for (let ii = 0; ii < 4; ii++) {
        patch.xplus[ii].xghost = patch.xplus[ii].x + tt*buffer;
        patch.xplus[ii].yghost = patch.xplus[ii].y;
      }
      for (let particle of patch.xplusparticles) {
        particle.xghostxp = particle.xnew + tt*buffer;
        particle.yghostxp = particle.ynew;
      }
    }

    function moveYMinus(patch, tt) {
      for (let ii = 0; ii < 4; ii++) {
        patch.yminus[ii].xghost = patch.yminus[ii].x;
        patch.yminus[ii].yghost = patch.yminus[ii].y - tt*buffer;
      }
      for (let particle of patch.yminusparticles) {
        particle.xghostym = particle.xnew;
        particle.yghostym = particle.ynew - tt*buffer;
      }
      for (let particle of patch.yminusxminusparticles) {
        particle.xghostymxm = particle.xghostxm;
        particle.yghostymxm = particle.yghostxm - tt*buffer;
      }
      for (let particle of patch.yminusxplusparticles) {
        particle.xghostymxp = particle.xghostxp;
        particle.yghostymxp = particle.yghostxp - tt*buffer;
      }
    }

    function moveYPlus(patch, tt) {
      for (let ii = 0; ii < 4; ii++) {
        patch.yplus[ii].xghost = patch.yplus[ii].x;
        patch.yplus[ii].yghost = patch.yplus[ii].y + tt*buffer;
      }
      for (let particle of patch.yplusparticles) {
        particle.xghostyp = particle.xnew;
        particle.yghostyp = particle.ynew + tt*buffer;
      }
      for (let particle of patch.yplusxminusparticles) {
        particle.xghostypxm = particle.xghostxm;
        particle.yghostypxm = particle.yghostxm + tt*buffer;
      }
      for (let particle of patch.yplusxplusparticles) {
        particle.xghostypxp = particle.xghostxp;
        particle.yghostypxp = particle.yghostxp + tt*buffer;
      }
    }

    function createUpdatedParticleLists() {
      for (let patch of patches) {
        patch.totalparticles = patch.totalparticles.concat(patch.particles);
        for (let particle of particles) {
          particle.xinit = particle.xnew;
          particle.yinit = particle.ynew;
          let pp = {x: particle.xghostxm, y:particle.yghostxm};
          if (inside(pp, patchesExtended[patch.id])) {
            let partcopy = Object.assign({}, particle);
            partcopy.xinit = particle.xghostxm;
            partcopy.yinit = particle.yghostxm;
            patch.totalparticles.push(partcopy);
          }
          pp = {x: particle.xghostxp, y:particle.yghostxp};
          if (inside(pp, patchesExtended[patch.id])) {
            let partcopy = Object.assign({}, particle);
            partcopy.xinit = particle.xghostxp;
            partcopy.yinit = particle.yghostxp;
            patch.totalparticles.push(partcopy);
          }
          pp = {x: particle.xghostym, y:particle.yghostym};
          if (inside(pp, patchesExtended[patch.id])) {
            let partcopy = Object.assign({}, particle);
            partcopy.xinit = particle.xghostym;
            partcopy.yinit = particle.yghostym;
            patch.totalparticles.push(partcopy);
          }
          pp = {x: particle.xghostyp, y:particle.yghostyp};
          if (inside(pp, patchesExtended[patch.id])) {
            let partcopy = Object.assign({}, particle);
            partcopy.xinit = particle.xghostyp;
            partcopy.yinit = particle.yghostyp;
            patch.totalparticles.push(partcopy);
          }
          pp = {x: particle.xghostymxm, y:particle.yghostymxm};
          if (inside(pp, patchesExtended[patch.id])) {
            let partcopy = Object.assign({}, particle);
            partcopy.xinit = particle.xghostymxm;
            partcopy.yinit = particle.yghostymxm;
            patch.totalparticles.push(partcopy);
          }
          pp = {x: particle.xghostymxp, y:particle.yghostymxp};
          if (inside(pp, patchesExtended[patch.id])) {
            let partcopy = Object.assign({}, particle);
            partcopy.xinit = particle.xghostymxp;
            partcopy.yinit = particle.yghostymxp;
            patch.totalparticles.push(partcopy);
          }
          pp = {x: particle.xghostypxm, y:particle.yghostypxm};
          if (inside(pp, patchesExtended[patch.id])) {
            let partcopy = Object.assign({}, particle);
            partcopy.xinit = particle.xghostypxm;
            partcopy.yinit = particle.yghostypxm;
            patch.totalparticles.push(partcopy);
          }
          pp = {x: particle.xghostypxp, y:particle.yghostypxp};
          if (inside(pp, patchesExtended[patch.id])) {
            let partcopy = Object.assign({}, particle);
            partcopy.xinit = particle.xghostypxp;
            partcopy.yinit = particle.yghostypxp;
            patch.totalparticles.push(partcopy);
          }
        }
      }
    }


    // Simulate motion of particles in each patch
    var simulateMotion = (tt) => {
      for (let patch of patches) {
        for (let part of patch.totalparticles) {
          part.xnew = part.xinit + part.xvel*tt;
          part.ynew = part.yinit + part.yvel*tt;
        }
      }
    }

    function doAnimation() {

      
      const maxTime = 3000;
      const ease = d3.easeCubic;
      let simulationX = d3.timer(
        function(elapsed) {
          const t = Math.min(1, ease(elapsed/maxTime));
   
          simulateMotion(t);

          context.clearRect(0, 0, width, height);
          drawPatches();
          drawMovedGhostRegionsX();
          drawMovedGhostRegionsY();
          drawDisplacedParticles();

          if (t === 1) {
            simulationX.stop();
          }
        });

    }

    function reset() {
      particles = [];
      patches = [];
      patchesExtended = [];
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

      for (let part of particles) {
        context.beginPath()
        context.arc(xscale(part.xnew), yscale(part.ynew), 5, 0, 2.0*Math.PI); 
        //context.fillStyle = "#aa0902";
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();
      }
    }

    function drawDisplacedParticles() {

      let patchID = 0;
      for (let patch of patches) {
        for (let part of patch.totalparticles) {
          context.beginPath()
          if (patch.id != part.patch) {
            context.arc(xscale(part.xnew), yscale(part.ynew), 5, 0, 1.5*Math.PI, false); 
            context.lineTo(xscale(part.xnew), yscale(part.ynew));
            context.lineTo(xscale(part.xnew)+5, yscale(part.ynew));
          } else {
            context.arc(xscale(part.xnew), yscale(part.ynew), 5, 0, 2.0*Math.PI); 
          }
          //context.fillStyle = "#aa0902";
          context.fillStyle = patchColors[part.patch];
          //context.fillStyle = patchColors[patchID];
          context.fill();
          context.stroke();
        }
        patchID++;
      }
    }

    function drawMovedGhostRegionsX() {
      context.save();
      context.lineWidth = 2;
      context.globalAlpha = 0.2;
      let patchID = 0;
      for (let patch of patches) {
        context.fillStyle = patchColors[patchID];
        if (patch[4].x != "xminus") {
          drawBox(patch.xminus, true, "ghost");
        }
        if (patch[4].x != "xplus") {
          drawBox(patch.xplus, true, "ghost");
        }
        patchID++;
      }
      context.restore();
    }

    function drawMovedGhostRegionsY(init = false) {
      let ghost = "ghost";
      if (init) {
        ghost = "none";
      }
      context.save();
      context.lineWidth = 2;
      context.globalAlpha = 0.2;
      let patchID = 0;
      for (let patch of patches) {
        context.fillStyle = patchColors[patchID];
        if (patch[4].y != "yminus") {
          drawBox(patch.yminus, true, ghost);
        }
        if (patch[4].y != "yplus") {
          drawBox(patch.yplus, true, ghost);
        }
        patchID++;
      }
      context.restore();
    }

    function drawMovedGhostParticlesX() {

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
      }
    }

    function drawCopiedGhostParticlesX() {

      for (let part of particles) {
        context.beginPath()
        context.arc(xscale(part.xghostymxm), yscale(part.yghostymxm), 5, 0, 2.0*Math.PI); 
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();

        context.beginPath()
        context.arc(xscale(part.xghostymxp), yscale(part.yghostymxp), 5, 0, 2.0*Math.PI); 
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();

        context.beginPath()
        context.arc(xscale(part.xghostypxm), yscale(part.yghostypxm), 5, 0, 2.0*Math.PI); 
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();

        context.beginPath()
        context.arc(xscale(part.xghostypxp), yscale(part.yghostypxp), 5, 0, 2.0*Math.PI); 
        context.fillStyle = patchColors[part.patch];
        context.fill();
        context.stroke();
      }
    }

    function drawMovedGhostParticlesY() {

      for (let part of particles) {
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

    let run = function () {
        initialize();
        doAnimation();
    }
    let restartAnimation =  function () {
        reset();
        initialize();
        doAnimation();
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

  particleMotion.run();


  /*
  function drawHiddenCanvas() {
    drawMovedParticles();
    drawMovedGhostRegionsX();
    drawMovedGhostParticlesX();
    drawCopiedGhostParticlesX(1.0);
    drawMovedGhostRegionsY();
    drawMovedGhostParticlesY();
  }

  function drawTooltip() {
    d3.select("#particle-migrate-plimpton") .on('mousemove', function() {
      drawHiddenCanvas();
      var mouseX = d3.event.layerX || d3.event.offsetX;
      var mouseY = d3.event.layerY || d3.event.offsetY;
      var hiddenCtx = hiddenCanvas.node().getContext('2d');
      var xx = xscale.invert(mouseX);
      var yy = yscale.invert(mouseY);
      var str = "x: " + xx + " y: " + yy;
      d3.select('#tooltip')
        .style('opacity', 0.8)
        .style('top', d3.event.pageY + 5 + 'px')
        .style('left', d3.event.pageX + 5 + 'px')
        .html(str);
    });
  }
  */

  /*
  run();
  drawTooltip();
  */
    
  


