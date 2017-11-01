function testCircleLine

  sides();
%  corners();

end

function sides()
  figure
  xplus();
  %xminus();
  %yminus();
end

function corners()
  figure
  xplusyminus();
  %xplusyplus();
  %xminusyplus();
  %xminusyminus();
end

function xplus()
  angle = 90*pi/180;
  cx = 2; cy = 1.5; lx = 1; ly = 2;
  mx = 0; my = 0; ex = 1; ey = 3;
  [u0, u1, N] = drawEllipse(cx, cy, lx, ly, angle);
  [v0, v1] = drawRect(mx, my, cx, cy, ex, ey, 0.5*angle, N);
end

function xminus()
  angle = 0*pi/180;
  cx = -2; cy = -1.5; lx = 1; ly = 2;
  mx = 0; my = 0; ex = 1; ey = 3;
  [u0, u1, N] = drawEllipse(cx, cy, lx, ly, angle);
  [v0, v1] = drawRect(mx, my, cx, cy, ex, ey, angle, N);
end

function yminus()
  angle = 90*pi/180;
  cx = -0.5; cy = -2.5; lx = 1; ly = 2;
  mx = 0; my = 0; ex = 1; ey = 3;
  [u0, u1, N] = drawEllipse(cx, cy, lx, ly, angle);
  [v0, v1] = drawRect(mx, my, cx, cy, ex, ey, angle, N);
end

function xplusyminus()
  angle = -40*pi/180;
  cx = 2; cy = -2; lx = 1; ly = 2;
  mx = 0; my = 0; ex = 1; ey = 3;
  [u0, u1, N] = drawEllipse(cx, cy, lx, ly, angle);
  [v0, v1] = drawRect(mx, my, cx, cy, ex, ey, 1.5*angle, N);
end

function xplusyplus()
  angle = 40*pi/180;
  cx = 2; cy = 2; lx = 1; ly = 2;
  mx = 0; my = 0; ex = 1; ey = 3;
  [u0, u1, N] = drawEllipse(cx, cy, lx, ly, angle);
  [v0, v1] = drawRect(mx, my, cx, cy, ex, ey, angle, N);
end

function xminusyplus()
  angle = -20*pi/180;
  cx = -2; cy = 2; lx = 1; ly = 2;
  mx = 0; my = 0; ex = 1; ey = 3;
  [u0, u1, N] = drawEllipse(cx, cy, lx, ly, angle);
  [v0, v1] = drawRect(mx, my, cx, cy, ex, ey, angle, N);
end

function xminusyminus()
  angle = 20*pi/180;
  cx = -2; cy = -2; lx = 1; ly = 2;
  mx = 0; my = 0; ex = 1; ey = 3;
  [u0, u1, N] = drawEllipse(cx, cy, lx, ly, angle);
  [v0, v1] = drawRect(mx, my, cx, cy, ex, ey, angle, N);
end

function [u0, u1, N] = drawEllipse(cx, cy, ax, ay, angle)

  theta = linspace(0, 360*pi/180, 30);
  x = ax*cos(theta);
  y = ay*sin(theta);
  M = [cos(angle) -sin(angle);sin(angle) cos(angle)];
  xx =  [x;y];
  xxp = M*xx;
  plot(xxp(1,:) + cx, xxp(2,:) + cy, 'r--', 'LineWidth', 1); hold on;
  plot(cx, cy, 'r.', 'MarkerSize', 10); hold on;
  axis equal;

  u0o = [0 ax; 0 0];
  u1o = [0 0; 0 ay];
  u0p = M*u0o;
  u1p = M*u1o;
  plot(u0p(1,:) + cx, u0p(2,:) + cy, 'r--', 'LineWidth', 1); hold on;
  plot(u1p(1,:) + cx, u1p(2,:) + cy, 'r--', 'LineWidth', 1); hold on;

  u0 = (u0p(:,2) - u0p(:,1));
  u1 = (u1p(:,2) - u1p(:,1));
  u0 = u0/norm(u0);
  u1 = u1/norm(u1);

  N = [u0(1,1)/ax u0(2,1)/ax; u1(1,1)/ay u1(2,1)/ay];
  xxpp = N*xxp;
  plot(xxpp(1,:), xxpp(2,:), 'r-', 'LineWidth', 2); hold on;

end
  
function [v0s, v1s] = drawRect(mx, my, cx, cy, ax, ay, angle, N)

  x = [mx-ax mx+ax mx+ax mx-ax mx-ax];
  y = [my-ay my-ay my+ay my+ay my-ay];
  M = [cos(angle) -sin(angle);sin(angle) cos(angle)];
  xx =  [x;y];
  xxp = M*xx;

  plot(xxp(1,:) + mx, xxp(2,:) + my, 'b--', 'LineWidth', 1); hold on;
  plot(mx, my, 'b.', 'MarkerSize', 10); hold on;

  v0o = [0 ax; 0 0];
  v1o = [0 0; 0 ay];
  v0p = M*v0o;
  v1p = M*v1o;
  plot(v0p(1,:) + mx, v0p(2,:) + my, 'b--', 'LineWidth', 1); hold on;
  plot(v1p(1,:) + mx, v1p(2,:) + my, 'g--', 'LineWidth', 1); hold on;

  v0 = (v0p(:,2) - v0p(:,1));
  v1 = (v1p(:,2) - v1p(:,1));
  v0 = v0/norm(v0);
  v1 = v1/norm(v1);

  ww = [mx - cx; mx - cy];
  mp = N*ww
  v0s = N*(v0*ax);
  v1s = N*(v1*ay);
  v0c = v0s/norm(v0s) + mp;
  v1c = v1s/norm(v1s) + mp;

  plot(mp(1,1), mp(2,1), 'b.', 'MarkerSize', 10); hold on;
  plot([mp(1,1) v0c(1,1)], [mp(2,1) v0c(2,1)], 'b-', 'LineWidth', 2); hold on;
  plot([mp(1,1) v1c(1,1)], [mp(2,1) v1c(2,1)], 'g-', 'LineWidth', 2); hold on;

  xx = [mp - v0s*ax - v1s*ay ...
        mp + v0s*ax - v1s*ay ...
        mp + v0s*ax + v1s*ay ...
        mp - v0s*ax + v1s*ay ...
        mp - v0s*ax - v1s*ay];
  plot(xx(1,:), xx(2,:), 'b-', 'LineWidth', 2); hold on;

  [t1, d1] = distance(0, 0, xx(1,1:2), xx(2,1:2))
  [t2, d2] = distance(0, 0, xx(1,2:3), xx(2,2:3))
  [t3, d3] = distance(0, 0, xx(1,3:4), xx(2,3:4))
  [t4, d4] = distance(0, 0, xx(1,4:5), xx(2,4:5))

  p1 = xx(:,1) + (xx(:,2)-xx(:,1))*t1;
  plot([xx(1,1) p1(1,1)],[xx(2,1) p1(2,1)],'m-', 'LineWidth', 2);
  plot([0 p1(1,1)],[0 p1(2,1)],'m--', 'LineWidth', 2);

  p2 = xx(:,2) + (xx(:,3)-xx(:,2))*t2;
  plot([xx(1,2) p2(1,1)],[xx(2,2) p2(2,1)],'m-', 'LineWidth', 2);
  plot([0 p2(1,1)],[0 p2(2,1)],'m--', 'LineWidth', 2);

  p3 = xx(:,3) + (xx(:,4)-xx(:,3))*t3;
  plot([xx(1,3) p3(1,1)],[xx(2,3) p3(2,1)],'m-', 'LineWidth', 2);
  plot([0 p3(1,1)],[0 p3(2,1)],'m--', 'LineWidth', 2);

  p4 = xx(:,4) + (xx(:,5)-xx(:,4))*t4;
  plot([xx(1,4) p4(1,1)],[xx(2,4) p4(2,1)],'m-', 'LineWidth', 2);
  plot([0 p4(1,1)],[0 p4(2,1)],'m--', 'LineWidth', 2);
end

function [t, d] = distance(px, py, segx, segy)

  x0 = [px; py; 0];
  x1 = [segx(1); segy(1); 0];
  x2 = [segx(2); segy(2); 0];

  x1x0 = x1 - x0;
  x2x0 = x2 - x0;
  x2x1 = x2 - x1;
  

  t = -dot(x1x0, x2x1)/norm(x2x1)^2;
  d = norm(cross(-x1x0, -x2x0))/norm(x2x1);
end
