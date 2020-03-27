function test
  plotBench()
  for i=1:5
    plotJoint()
  end
end

function plotBench()
  plotBox([0.501 5],[15 15])
  plotBox([0.501 0.3],[30 5])
  plotBox([10.0 4.0],[10.2 15])
  
end

function plotBox(lo, hi)
  x = [lo(1) hi(1) hi(1) lo(1) lo(1)];
  y = [lo(2) lo(2) hi(2) hi(2) lo(2)];
  plot(x, y); hold on;
  axis equal;
end

function plotJoint()
  rad = 10;
  [x, y] = ginput(2);
  [x(1) y(1); x(2) y(2)]
  xmin = min(x);
  ymin = min(y);
  plot(x, y, 'x'); 
  [xp, yp, R] = rotateToX(x-xmin,y-ymin);
  [lo, hi] = plotJointY(xp, yp, rad);
  [p1, p2, p3, p4] = rotateBack(lo, hi, R);
  p1s = translate(p1, xmin, ymin);
  p2s = translate(p2, xmin, ymin);
  p3s = translate(p3, xmin, ymin);
  p4s = translate(p4, xmin, ymin);
  plot([p1s(1) p2s(1) p3s(1) p4s(1) p1s(1)],[p1s(2) p2s(2) p3s(2) p4s(2) p1s(2)]);
end

function [xp, yp, R] = rotateToX(x, y)

  %plot(x, y, 'o'); 
  
  theta = atan2(y(2)-y(1), x(2)-x(1));
  R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
  p1 = [x(1); y(1)];
  p2 = [x(2); y(2)];
  p1p = R'*p1;
  p2p = R'*p2;
  xp = [p1p(1) p2p(1)];
  yp = [p1p(2) p2p(2)];
  %plot(xp, yp, 's');
end

function [lo, hi] = plotJointY(xp, yp, rad)

  lo = [xp(1); yp(1) - rad];
  hi = [xp(2); yp(2) + rad];
  %plotBox(lo, hi);
end

function [p1p, p2p, p3p, p4p] = rotateBack(lo, hi, R)

  p1 = lo;
  p2 = [hi(1);lo(2)];
  p3 = hi;
  p4 = [lo(1);hi(2)];
  p1p = R*p1;
  p2p = R*p2;
  p3p = R*p3;
  p4p = R*p4;
  %plot([p1p(1) p2p(1) p3p(1) p4p(1) p1p(1)],[p1p(2) p2p(2) p3p(2) p4p(2) p1p(2)]);
end

function [ps] = translate(p, xmin, ymin)
  ps(1) = p(1) + xmin;
  ps(2) = p(2) + ymin;
end
