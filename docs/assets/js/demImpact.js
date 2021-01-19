  var canvasBB = document.querySelector("#ballball");
  var contextBB = canvasBB.getContext("2d");
  var widthBB = canvasBB.width;
  var heightBB = canvasBB.height;
  var theta = 2 * Math.PI;

  var balls = [{x: 100, y: heightBB/2, r: 50, vx: 1.0, vy: 0}, 
               {x: 300, y: heightBB/2, r: 30, vx: 0.0, vy: 0}];


  var simulation = d3.forceSimulation(balls)
    .velocityDecay(0.002)
    .force("collide", d3.forceCollide()
                        .radius(function(d) {return d.r + 0.5})
                        .iterations(1))
    .on("tick", drawBalls)

  d3.interval(function() {restart();}, 7000, d3.now());

  function restart() {
    balls = [{x: 100, y: heightBB/2, r: 50, vx: 1.0, vy: 0}, 
               {x: 300, y: heightBB/2, r: 30, vx: 0.0, vy: 0}];
    simulation.nodes(balls);
    simulation.alpha(1).restart();
  }

  function drawBalls() {
    contextBB.clearRect(0, 0, widthBB, heightBB);
    contextBB.save();
    contextBB.fillStyle = "#ddd";
    contextBB.strokeStyle = "#333";
    contextBB.beginPath();
    contextBB.arc(balls[0].x, balls[0].y, balls[0].r, 0, theta);
    contextBB.fill();
    contextBB.stroke();
    contextBB.beginPath();
    contextBB.arc(balls[1].x, balls[1].y, balls[1].r, 0, theta);
    contextBB.fill();
    contextBB.stroke();
    contextBB.restore();

    drawWall();
    drawSpring();
  }

  function drawWall() {
    let wallWidth = 50;
    let wallHeight = 200;
    let halfHeight = wallHeight/2;
    let xStart = 450;
    let yStart = heightBB/2 - halfHeight;
    contextBB.save();
    contextBB.translate(xStart, yStart);
    contextBB.fillStyle = "#ddd";
    contextBB.fillRect(0, 0, wallWidth, wallHeight);
    contextBB.restore();
  }

  function drawSpring() {

    let xStartInit = 300 - balls[1].r;
    let xStart0 = balls[0].x + balls[0].r;
    let xStart = balls[1].x-balls[1].r;
    let xEnd = 450;
    let yStart = heightBB/2;
    let yEnd = heightBB/2;

    contextBB.save();
    contextBB.translate(xStartInit, yStart);
    contextBB.fillRect(0, 70,2,30);
    contextBB.font = "16px Arial";
    contextBB.fillText("overlap", 10, 85);

    contextBB.translate(xStart0-xStartInit, 0);
    contextBB.fillRect(0,-100,2,30);
    contextBB.fillText("gap", 10, -90);

    contextBB.translate(xStart-xStart0, 0);
    contextBB.fillRect(-2,-30,5,60);

    contextBB.fillRect(0,-100,2,30);
    contextBB.fillRect(0, 70,2,30);
    drawArrow(0, -85, 50, -85);

    let springLen = xEnd - xStart;
    let numLoops = 20;
    let loopLen = springLen/(numLoops+2);

    contextBB.beginPath();
    contextBB.moveTo(0,0);
    contextBB.lineTo(loopLen, 0);

    let loopWidth = 40;
    let halfWidth = loopWidth/2;

    for (let ii = 0; ii < numLoops; ii++) {
      let xControl1 = loopLen * (1 + ii);
      let yControl1 = halfWidth;
      let xControl2 = loopLen * (1 + ii + 0.5);
      let yControl2 = yControl1;
      let xControl3 = xControl2;
      let yControl3 = 0;
      contextBB.bezierCurveTo(xControl1, yControl1, 
        xControl2, yControl2, xControl3, yControl3);
      xControl1 = xControl3;
      yControl1 = - halfWidth;
      xControl2 = loopLen * (1 + ii + 1.0);
      yControl2 = yControl1;
      xControl3 = xControl2;
      yControl3 = 0;
      contextBB.bezierCurveTo(xControl1, yControl1, 
        xControl2, yControl2, xControl3, yControl3);
    }

    contextBB.lineTo(springLen, 0);
    contextBB.stroke();
    contextBB.restore();
  }

  function drawArrow(xStart, yStart, xEnd, yEnd) {
    let headLen = 10.0;
    let angle = Math.PI/13.0;
    contextBB.moveTo(xStart, yStart);
    contextBB.lineTo(xEnd, yEnd);
    contextBB.moveTo(xEnd, yEnd);
    contextBB.lineTo(xEnd - headLen*Math.cos(-angle), yEnd - headLen*Math.sin(-angle));
    contextBB.moveTo(xEnd, yEnd);
    contextBB.lineTo(xEnd - headLen*Math.cos(angle), yEnd - headLen*Math.sin(angle));
    contextBB.fillText("d", xEnd+10, yEnd);
    contextBB.stroke()
  }
