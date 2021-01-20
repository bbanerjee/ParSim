  var canvas = document.querySelector("#particle-scatter");
  var context = canvas.getContext("2d");
  var width = canvas.width;
  var height = canvas.height;

  var xscale = d3.scaleLinear()
                 .domain([-1, 9.6])
                 .range([0, width]);
  var yscale = d3.scaleLinear()
                 .domain([-1, 9.6])
                 .range([height, 0]);

  var line = d3.line()
               .x(function(d) {return xscale(d.x);})
               .y(function(d) {return yscale(d.y);});

  var initSquare = [{x: 0, y: 0},
                    {x: 2, y: 0},
                    {x: 2, y: 2},
                    {x: 0, y: 2}];

  var particles = [];
  var patches = [];
  var patchMaps = [];
  var patchColors = colorbrewer.Set1[9];

  initialize();
  doAnimation();
  
  function initialize() {
    createParticles(15);
    createPatches();
    locateParticlesInPatch();
    moveParticlesToPatches();

    context.save();
    context.clearRect(0, 0, width, height);
    context.fillStyle = "#87ceeb";
    context.globalAlpha = 0.6;
    drawSquare(initSquare, true);
    let label = "P0";
    context.font = "16px Arial";
    context.fillStyle = "#aa0902";
    context.fillText(label, xscale(initSquare[0].x)+20, yscale(initSquare[0].y)-20);
    context.restore();

    drawAxes();

    setTimeout(drawParticles, 100);
    setTimeout(drawPatches, 100);
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

        if (t === 1) {
          simulation.stop();
        }
      });
  }

  function restartAnimation() {
    particles = [];
    patches = [];
    patchMaps = [];
    initialize();
    doAnimation();
  }


  function drawSquare(curSquare, dashed) {
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

  function drawAxes() {
    var xaxis = [{x: xscale(0), y: yscale(0)},
                 {x: xscale(3), y: yscale(0)}];
    var yaxis = [{x: xscale(0), y: yscale(0)},
                 {x: xscale(0), y: yscale(3)}];
    drawArrow(xaxis[0].x, xaxis[0].y, xaxis[1].x, xaxis[1].y, "x");
    drawArrow(yaxis[0].x, yaxis[0].y, yaxis[1].x, yaxis[1].y, "y");
  }

  function drawArrow(xStart, yStart, xEnd, yEnd, label) {
    let orient = Math.atan2(yEnd - yStart, xEnd - xStart);
    let headLen = 20.0;
    let angle = Math.PI/13.0;
    context.moveTo(xStart, yStart);
    context.lineTo(xEnd, yEnd);
    context.moveTo(xEnd, yEnd);
    context.lineTo(xEnd - headLen*Math.cos(orient-angle), 
                   yEnd - headLen*Math.sin(orient-angle));
    context.moveTo(xEnd, yEnd);
    context.lineTo(xEnd - headLen*Math.cos(orient+angle), 
                   yEnd - headLen*Math.sin(orient+angle));
    context.font = "16px Arial";
    context.fillStyle = "#aa0902";
    context.fillText(label, xEnd-20, yEnd+20);
    context.stroke();
  }

  // Create particles
  function createParticles(num) {
    let count = 0;
    while (count < num) {
      let xCen = d3.randomUniform(0.1, 1.9)();
      let yCen = d3.randomUniform(0.1, 1.9)();
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

  // Draw particles
  function drawParticles() {
    for (let part of particles) {
      if (!part.moved) {
        context.beginPath()
        context.arc(xscale(part.x), yscale(part.y), 5, 0, 2.0*Math.PI); 
        context.fillStyle = "#aa0902";
        context.fill();
        context.stroke();
      }
    }
  }

  // Create patches
  function createPatches() {

    let xShift = 3;
    let yShift = 3;
    let numPatch = 3;
    let patchWidth = 2; 
    let buffer = 0.2;
    let patchID = 0;
    for (let ii = 0; ii < numPatch; ii++) {
      let xmap = ii*2/3;
      let xloc = xShift + ii*(patchWidth+buffer);
      for (let jj = 0; jj < numPatch; jj++) {
        let yloc = yShift + jj*(patchWidth+buffer);
        let ymap = jj*2/3;
        patches.push([{x: xloc, y: yloc},
                      {x: xloc + patchWidth, y: yloc},
                      {x: xloc + patchWidth, y: yloc + patchWidth},
                      {x: xloc, y: yloc + patchWidth}]);
        patchMaps.push([{x: xmap, y: ymap},
                        {x: xmap + 2/3, y: ymap},
                        {x: xmap + 2/3, y: ymap + 2/3},
                        {x: xmap, y: ymap + 2/3}]);
      }
    }
  }

  // Draw patches
  function drawPatches() {
    let patchID = 0;
    for (let patch of patches) {
      context.save();
      context.lineWidth = 2;
      if (patchID === 0) {
        context.fillStyle = "#87ceeb";
      } else {
        context.fillStyle = "#9acd32";
      }
      context.globalAlpha = 0.6;
      drawSquare(patch, false);

      let label = "P" + patchID;
      context.font = "16px Arial";
      context.fillStyle = "#aa0902";
      context.fillText(label, xscale(patch[0].x)+20, yscale(patch[0].y)-20);
      context.stroke();
      context.restore();
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
        }
        patchID++;
      }
    }
  }

  // Draw moved particles
  function drawMovedParticles() {

    context.clearRect(0, 0, width, height);

    // Initial status
    context.save();
    context.fillStyle = "#87ceeb";
    context.strokeStyle = "#333";
    context.lineWidth = 2;
    context.globalAlpha = 0.6;
    drawSquare(initSquare, false);

    let label = "P0";
    let xpos = 0.5*(initSquare[0].x+initSquare[1].x);
    let ypos = initSquare[3].y;
    context.font = "16px Arial";
    context.fillStyle = "#aa0902";
    context.fillText(label, xscale(xpos), yscale(ypos)-5);
    context.restore();

    drawParticles();
    drawPatches();
    drawAxes();

    for (let part of particles) {
      context.beginPath()
      context.arc(xscale(part.xpos), yscale(part.ypos), 5, 0, 2.0*Math.PI); 
      //context.fillStyle = "#aa0902";
      context.fillStyle = patchColors[part.patch];
      context.fill();
      context.stroke();
    }
  }

