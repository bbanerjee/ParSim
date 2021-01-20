  var canvas = document.querySelector("#rot-elem");
  var context = canvas.getContext("2d");
  var width = canvas.width;
  var height = canvas.height;

  var xscale = d3.scaleLinear()
                 .domain([-4, 2])
                 .range([0, width]);
  var yscale = d3.scaleLinear()
                 .domain([-4, 2])
                 .range([height, 0]);

  var line = d3.line()
               .x(function(d) {return xscale(d.x);})
               .y(function(d) {return yscale(d.y);});

  var initSquare = [{x: -1, y: -1},
                    {x: 1, y: -1},
                    {x: 1, y: 1},
                    {x: -1, y: 1}];
  var square = JSON.parse(JSON.stringify(initSquare));
  var strainedSquare = JSON.parse(JSON.stringify(initSquare));

  var curAngle = 0;

  doSimulation();
  
  function doSimulation() {
    var simulation = d3.interval(function(elapsed) {
       rotAngle = 1;
       curAngle += rotAngle;
       rotateSquare(rotAngle);
       computeStrains(curAngle+1);
       if (curAngle >= 360) {
         simulation.stop();
       }
    }, 60);

    context.clearRect(0, 0, width, height);
    drawSquare(initSquare, true);
    drawSquare(square, false);
    drawAxes();
  }

  function restartSimulation() {
    reset();
    doSimulation();
  }

  function reset() {
    square = JSON.parse(JSON.stringify(initSquare));
    curAngle = 0;
  }

  function computeStrains(angle) {
    let theta = angle*Math.PI/180;
    let costheta = Math.cos(theta);
    let sintheta = Math.sin(theta);
    let epsx = costheta - 1.0;
    let epsy = costheta - 1.0;
    let epsxp = epsx*costheta*costheta + epsy*sintheta*sintheta;
    let epsyp = epsx*sintheta*sintheta + epsy*costheta*costheta;
    let deltaxp = 2*epsxp;
    let deltayp = 2*epsyp;
    strainedSquare = JSON.parse(JSON.stringify(initSquare));
    strainedSquare[1].x += deltaxp;
    strainedSquare[2].x += deltaxp;
    strainedSquare[2].y += deltayp;
    strainedSquare[3].y += deltayp;
    let originx = strainedSquare[0].x;
    let originy = strainedSquare[0].y;
    for (var vertex of strainedSquare) {
      let transx = vertex.x - originx;
      let transy = vertex.y - originy;
      let rotx = transx*Math.cos(theta) + transy*Math.sin(theta);
      let roty = -transx*Math.sin(theta) + transy*Math.cos(theta);
      vertex.x = rotx + originx;
      vertex.y = roty + originy;
    } 
  }

  function rotateSquare(angle) {
    let theta = angle*Math.PI/180;
    let originx = square[0].x;
    let originy = square[0].y;
    for (var vertex of square) {
      let transx = vertex.x - originx;
      let transy = vertex.y - originy;
      let rotx = transx*Math.cos(theta) + transy*Math.sin(theta);
      let roty = -transx*Math.sin(theta) + transy*Math.cos(theta);
      vertex.x = rotx + originx;
      vertex.y = roty + originy;
    } 

    context.clearRect(0, 0, width, height);

    // Initial status
    context.save();
    context.fillStyle = "#87ceeb";
    context.strokeStyle = "#333";
    context.lineWidth = 1;
    context.globalAlpha = 0.6;
    drawSquare(initSquare, true);
    context.fillRect(20, 20, 20, 20);
    context.fillText("Initial", 50, 35);
    context.restore();

    // Strained status
    context.save();
    context.fillStyle = "#f4a460";
    context.lineWidth = 1;
    context.globalAlpha = 0.6;
    drawSquare(strainedSquare, false);
    context.fillRect(20, 50, 20, 20);
    context.fillText("Strained", 50, 65);
    context.restore();

    // Rotated status
    context.save();
    context.lineWidth = 2;
    context.fillStyle = "#9acd32";
    context.globalAlpha = 0.6;
    drawSquare(square, false);
    context.fillRect(20, 80, 20, 20);
    context.fillText("Rotated", 50, 95);
    context.restore();
    drawAxes();
    drawAngle(xscale(square[0].x), yscale(square[0].y), 100, curAngle);
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
    var xaxis = [{x: xscale(-4), y: yscale(-1.0)},
                 {x: xscale(2.0), y: yscale(-1.0)}];
    var yaxis = [{x: xscale(-1.0), y: yscale(-4)},
                 {x: xscale(-1.0), y: yscale(2.0)}];
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

  function drawAngle(xCen, yCen, radius, angle) {
    let theta = angle*Math.PI/180;
    context.beginPath()
    context.setLineDash([5,5]);
    context.arc(xCen, yCen, radius, 0, theta); 
    context.font = "16px Arial";
    context.fillStyle = "#aa0902";
    context.fillText("\u03b2", xCen+radius/2, yCen+20);
    context.stroke();
    context.setLineDash([0]);
  }


