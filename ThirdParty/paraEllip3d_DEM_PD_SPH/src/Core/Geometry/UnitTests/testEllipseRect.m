function testEllipseRect

  figure;
  oneVisible();

  figure;
  twoVisible();

end

function oneVisible()
  xplus();
  xminus();
  yminus();
end

function xplus()

  cx = 2; cy = 1.5; lx = 1; ly = 2;
  [u0, u1] = drawEllipse(cx, cy, lx, ly, 0*pi/180);

  mx = 0; my = 0; ex = 1; ey = 3;
  [v0, v1] = drawRect(mx, my, ex, ey);

  drawLines(cx, cy, mx, my, v0, v1); 

end

function xminus()

  cx = -2; cy = -1.5; lx = 1; ly = 2;
  [u0, u1] = drawEllipse(cx, cy, lx, ly, 0*pi/180);

  mx = 0; my = 0; ex = 1; ey = 3;
  [v0, v1] = drawRect(mx, my, ex, ey);

  drawLines(cx, cy, mx, my, v0, v1); 

end

function yminus()

  cx = -0.5; cy = -2.5; lx = 1; ly = 2;
  [u0, u1] = drawEllipse(cx, cy, lx, ly, 90*pi/180);

  mx = 0; my = 0; ex = 1; ey = 3;
  [v0, v1] = drawRect(mx, my, ex, ey);

  drawLines(cx, cy, mx, my, v0, v1); 

end

function twoVisible()
  xplusyminus();
  xplusyplus();
  xminusyplus();
  xminusyminus();
end

function xplusyminus()
  cx = 2; cy = -2; lx = 1; ly = 2;
  [u0, u1] = drawEllipse(cx, cy, lx, ly, -40*pi/180);

  mx = 0; my = 0; ex = 1; ey = 3;
  [v0, v1] = drawRect(mx, my, ex, ey);

  drawLines(cx, cy, mx, my, v0, v1); 

end

function xplusyplus()
  cx = 2; cy = 2; lx = 1; ly = 2;
  [u0, u1] = drawEllipse(cx, cy, lx, ly, 40*pi/180);

  mx = 0; my = 0; ex = 1; ey = 3;
  [v0, v1] = drawRect(mx, my, ex, ey);

  drawLines(cx, cy, mx, my, v0, v1); 

end

function xminusyplus()
  cx = -2; cy = 2; lx = 1; ly = 2;
  [u0, u1] = drawEllipse(cx, cy, lx, ly, -20*pi/180);

  mx = 0; my = 0; ex = 1; ey = 3;
  [v0, v1] = drawRect(mx, my, ex, ey);

  drawLines(cx, cy, mx, my, v0, v1); 

end

function xminusyminus()
  cx = -2; cy = -2; lx = 1; ly = 2;
  [u0, u1] = drawEllipse(cx, cy, lx, ly, 20*pi/180);

  mx = 0; my = 0; ex = 1; ey = 3;
  [v0, v1] = drawRect(mx, my, ex, ey);

  drawLines(cx, cy, mx, my, v0, v1); 

end

function [u0, u1] = drawEllipse(cx, cy, ax, ay, angle)

  theta = linspace(0, 360*pi/180, 30);
  x = ax*cos(theta);
  y = ay*sin(theta);
  M = [cos(angle) -sin(angle);sin(angle) cos(angle)];
  xx =  [x;y];
  xxp = M*xx;
  plot(xxp(1,:) + cx, xxp(2,:) + cy, 'r-', 'LineWidth', 2); hold on;
  plot(cx, cy, 'r.', 'MarkerSize', 10); hold on;
  axis equal;

  u0o = [0 ax; 0 0];
  u1o = [0 0; 0 ay];
  u0p = M*u0o;
  u1p = M*u1o;
  plot(u0p(1,:) + cx, u0p(2,:) + cy, 'r-', 'LineWidth', 2); hold on;
  plot(u1p(1,:) + cx, u1p(2,:) + cy, 'r-', 'LineWidth', 2); hold on;

  u0 = (u0p(:,2) - u0p(:,1));
  u1 = (u1p(:,2) - u1p(:,1));
  u0 = u0/norm(u0);
  u1 = u1/norm(u1);

end
  
function [v0, v1] = drawRect(cx, cy, ax, ay)

  x = [cx-ax cx+ax cx+ax cx-ax cx-ax];
  y = [cy-ay cy-ay cy+ay cy+ay cy-ay];
  plot(x, y, 'b-', 'LineWidth', 2); hold on;
  plot(cx, cy, 'b.', 'MarkerSize', 10); hold on;

  v0o = [0 ax; 0 0];
  v1o = [0 0; 0 ay];
  plot(v0o(1,:) + cx, v0o(2,:) + cy, 'b-', 'LineWidth', 2); hold on;
  plot(v1o(1,:) + cx, v1o(2,:) + cy, 'b-', 'LineWidth', 2); hold on;

  v0 = (v0o(:,2) - v0o(:,1));
  v1 = (v1o(:,2) - v1o(:,1));
  v0 = v0/norm(v0);
  v1 = v1/norm(v1);

end

function drawLines(cx, cy, mx, my, v0, v1) 

  w = [cx - mx; cy - my];
  plot([mx cx], [my cy], 'm-', 'LineWidth', 2);

  d0 = dot(v0, w);
  d1 = dot(v1, w);
  d0v0 = d0*v0;
  d1v1 = d1*v1;
  plot([mx d0v0(1)],[my d0v0(2)], 'g--', 'LineWidth', 2);
  plot([mx d1v1(1)],[my d1v1(2)], 'g--', 'LineWidth', 2);
  plot([cx d0v0(1)],[cy d0v0(2)], 'g--', 'LineWidth', 2);
  plot([cx d1v1(1)],[cy d1v1(2)], 'g--', 'LineWidth', 2);

  [d0 d1]

end
