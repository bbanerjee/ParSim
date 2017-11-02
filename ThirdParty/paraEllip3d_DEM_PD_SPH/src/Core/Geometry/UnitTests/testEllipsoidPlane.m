function testEllipsoidPlane

  graphics_toolkit("gnuplot");
  %sides();
  corners();


end

function sides()
  %figure
  %xplus();
  %axis('equal');
  %figure
  %xminus();
  %axis('equal');
  %figure
  %yminus();
  %axis('equal');
end

function corners()
  %figure
  xplusyminus();
  axis('equal');
  %figure
  %xplusyplus();
  %axis('equal');
  %figure
  %xminusyplus();
  %axis('equal');
  %figure
  %xminusyminus();
  %axis('equal');
end

function xplus()
  axis = [0 0 1];
  angle = 90*pi/180;
  cx = 2; cy = 1.5; cz = 1.0; 
  lx = 1; ly = 2; lz = 3;
  mx = 0; my = 0; mz = 0; 
  ex = 1; ey = 3; ez = 2;
  cc = [cx; cy; cz];
  mm = [mx; my; mz];
  ll = [lx; ly; lz];
  ee = [ex; ey; ez];
  fig1 = figure;
  fig2 = figure;
  [u0, u1, u2, N] = drawEllipsoid(cc, ll, axis, angle, fig1, fig2);
  [v0, v1, v2] = drawOrientedBox(mm, cc, ee, axis, 0.5*angle, N, fig1, fig2);
end

function xminus()
  axis = [0 0 1];
  angle = 0*pi/180;
  cx = -2; cy = -1.5; cz = 1; 
  lx = 1; ly = 2; lz = 3;
  mx = 0; my = 0; mz = 0; 
  ex = 1; ey = 3; ez = 2;
  cc = [cx; cy; cz];
  mm = [mx; my; mz];
  ll = [lx; ly; lz];
  ee = [ex; ey; ez];
  fig1 = figure;
  fig2 = figure;
  [u0, u1, u2, N] = drawEllipsoid(cc, ll, axis, angle, fig1, fig2);
  [v0, v1, v2] = drawOrientedBox(mm, cc, ee, axis, angle, N, fig1, fig2);
end

function yminus()
  axis = [0 0 1];
  angle = 90*pi/180;
  cx = -0.5; cy = -2.5; cz = 1; 
  lx = 1; ly = 2; lz = 3;
  mx = 0; my = 0; mz = 0;
  ex = 1; ey = 3; ez = 2;
  cc = [cx; cy; cz];
  mm = [mx; my; mz];
  ll = [lx; ly; lz];
  ee = [ex; ey; ez];
  fig1 = figure;
  fig2 = figure;
  [u0, u1, u2, N] = drawEllipsoid(cc, ll, axis, angle, fig1, fig2);
  [v0, v1, v2] = drawOrientedBox(mm, cc, ee, axis, angle, N, fig1, fig2);
end

function xplusyminus()
  axis = [0 0 1];
  angle = -40*pi/180;
  cx = 2.6; cy = -2; cz = 1;
  lx = 1; ly = 2; lz = 3;
  mx = 0; my = 0; mz = 0;
  ex = 1; ey = 3; ez = 2;
  cc = [cx; cy; cz];
  mm = [mx; my; mz];
  ll = [lx; ly; lz];
  ee = [ex; ey; ez];
  fig1 = figure;
  fig2 = figure;
  [u0, u1, u2, N] = drawEllipsoid(cc, ll, axis, angle, fig1, fig2);
  [v0, v1, v2] = drawOrientedBox(mm, cc, ee, axis, 0, N, fig1, fig2);
end

function xplusyplus()
  axis = [0 0 1];
  angle = 40*pi/180;
  cx = 2; cy = 2; cz = 1;
  lx = 1; ly = 2; lz = 3;
  mx = 0; my = 0; mz = 0;
  ex = 1; ey = 3; ez = 2;
  cc = [cx; cy; cz];
  mm = [mx; my; mz];
  ll = [lx; ly; lz];
  ee = [ex; ey; ez];
  fig1 = figure;
  fig2 = figure;
  [u0, u1, u2, N] = drawEllipsoid(cc, ll, axis, angle, fig1, fig2);
  [v0, v1, v2] = drawOrientedBox(mm, cc, ee, axis, 0, N, fig1, fig2);
end

function xminusyplus()
  axis = [0 0 1];
  angle = -20*pi/180;
  cx = -2; cy = 2; cz = 1;
  lx = 1; ly = 2; lz = 3;
  mx = 0; my = 0; mz = 0;
  ex = 1; ey = 3; ez = 2;
  cc = [cx; cy; cz];
  mm = [mx; my; mz];
  ll = [lx; ly; lz];
  ee = [ex; ey; ez];
  fig1 = figure;
  fig2 = figure;
  [u0, u1, u2, N] = drawEllipsoid(cc, ll, axis, angle, fig1, fig2);
  [v0, v1, v2] = drawOrientedBox(mm, cc, ee, axis, 0, N, fig1, fig2);
end

function xminusyminus()
  axis = [0 0 1];
  angle = 20*pi/180;
  cx = -2; cy = -2; cz = 1;
  lx = 1; ly = 2; lz = 3;
  mx = 0; my = 0; mz = 0;
  ex = 1; ey = 3; ez = 2;
  cc = [cx; cy; cz];
  mm = [mx; my; mz];
  ll = [lx; ly; lz];
  ee = [ex; ey; ez];
  fig1 = figure;
  fig2 = figure;
  [u0, u1, u2, N] = drawEllipsoid(cc, ll, axis, angle, fig1, fig2);
  [v0, v1, v2] = drawOrientedBox(mm, cc, ee, axis, 0, N, fig1, fig2);
end

function [u0, u1, u2, N] = drawEllipsoid(cc, ll, axis, angle, fig1, fig2)

  k = 5;
  n = 2^k-1;
  theta = pi*(-n:2:n)/n;
  phi = (pi/2)*(-n:2:n)'/n;
  X = ll(1)*cos(phi)*cos(theta);
  Y = ll(2)*cos(phi)*sin(theta);
  Z = ll(3)*sin(phi)*ones(size(theta));
  [ntheta, nphi] = size(X)
  X = reshape(X, 1, nphi*ntheta); 
  Y = reshape(Y, 1, nphi*ntheta); 
  Z = reshape(Z, 1, nphi*ntheta); 
  xx = [X; Y; Z];

  M = [cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0; 0 0 1];
  xp = M*xx + cc;
  X = reshape(xp(1,:), ntheta, nphi);
  Y = reshape(xp(2,:), ntheta, nphi);
  Z = reshape(xp(3,:), ntheta, nphi);
  figure(fig1)
  hh = surf(X, Y, Z, 'edgecolor', 'none', 'facecolor', 'r', 'facealpha', 0.1); hold on;
  set (gca (), "plotboxaspectratio", [1 1 1])
  view(2)
  xlabel('x');
  ylabel('y');
  %axis('equal');

  u0o = [0 ll(1); 0 0; 0 0];
  u1o = [0 0; 0 ll(2); 0 0];
  u2o = [0 0; 0 0; 0 ll(3)];
  u0p = M*u0o + cc;
  u1p = M*u1o + cc;
  u2p = M*u2o + cc;
  plot3(u0p(1,:), u0p(2,:), u0p(3,:), 'r--', 'LineWidth', 1); hold on;
  plot3(u1p(1,:), u1p(2,:), u1p(3,:), 'r--', 'LineWidth', 1); hold on;
  plot3(u2p(1,:), u2p(2,:), u2p(3,:), 'r--', 'LineWidth', 1); hold on;
  plot3(cc(1), cc(2), cc(3), 'r.', 'MarkerSize', 12);

  u0 = (u0p(:,2) - u0p(:,1));
  u1 = (u1p(:,2) - u1p(:,1));
  u2 = (u2p(:,2) - u2p(:,1));
  u0 = u0/norm(u0);
  u1 = u1/norm(u1);
  u2 = u2/norm(u2);
  [u0'; u1'; u2']
  ppp = createBoxVertices(cc, ll, u0, u1, u2);
  drawBoxFaces(ppp, 'b');

  N = [u0(1,1)/ll(1) u0(2,1)/ll(1) u0(3,1)/ll(1); ...
       u1(1,1)/ll(2) u1(2,1)/ll(2) u1(3,1)/ll(2); ...
       u2(1,1)/ll(3) u2(2,1)/ll(3) u2(3,1)/ll(3)];

  X = reshape(X, 1, nphi*ntheta); 
  Y = reshape(Y, 1, nphi*ntheta); 
  Z = reshape(Z, 1, nphi*ntheta); 
  xx = [X; Y; Z];
  xxp = N*(xx - cc);
  X = reshape(xxp(1,:), ntheta, nphi);
  Y = reshape(xxp(2,:), ntheta, nphi);
  Z = reshape(xxp(3,:), ntheta, nphi);
  figure(fig2)
  hh = surf(X, Y, Z, 'edgecolor', 'none', 'facealpha', 0.5); hold on;
  plot3(0, 0, 0, 'r.', 'MarkerSize', 12);
  set (gca (), "plotboxaspectratio", [1 1 1])
  view(2)
  xlabel('x');
  ylabel('y');

end
  
function [v0s, v1s, v2s] = drawOrientedBox(mm, cc, ee, axis, angle, N, fig1, fig2)

  ax1 = [1;0;0]; ax2 = [0;1;0]; ax3 = [0;0;1];
  p = createBoxVertices(mm, ee, ax1, ax2, ax3);

  M = [cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0; 0 0 1];
  for ii=1:8
    pp(ii).xx = M*p(ii).xx;
  end

  figure(fig1);
  drawBoxFaces(pp, 'r');

  plot3(mm(1,1), mm(2,1), mm(3,1), 'b.', 'MarkerSize', 10); hold on;

  v0o = [0 ee(1); 0 0; 0 0];
  v1o = [0 0; 0 ee(2); 0 0];
  v2o = [0 0; 0 0; 0 ee(3)];
  v0p = M*v0o;
  v1p = M*v1o;
  v2p = M*v2o;
  plot3(v0p(1,:) + mm(1), v0p(2,:) + mm(2), v0p(3,:) + mm(3), 'b--', 'LineWidth', 1); hold on;
  plot3(v1p(1,:) + mm(1), v1p(2,:) + mm(2), v1p(3,:) + mm(3), 'b--', 'LineWidth', 1); hold on;
  plot3(v2p(1,:) + mm(1), v2p(2,:) + mm(2), v2p(3,:) + mm(3), 'b--', 'LineWidth', 1); hold on;

  v0 = (v0p(:,2) - v0p(:,1));
  v1 = (v1p(:,2) - v1p(:,1));
  v2 = (v2p(:,2) - v2p(:,1));
  v0 = v0/norm(v0);
  v1 = v1/norm(v1);
  v2 = v2/norm(v2);
  [v0'; v1';  v2']

  ww = mm - cc;
  mp = N*ww;
  v0s = N*(v0*ee(1));
  v1s = N*(v1*ee(2));
  v2s = N*(v2*ee(3));

  figure(fig2);
  %for i=1:8
  % ppp(i).xx = N*(pp(i).xx - cc);
  %end
  ee = [1;1;1];
  ppp = createBoxVertices(mp, ee, v0s, v1s, v2s);
  drawBoxFaces(ppp, 'm');

  %for i=1:length(ppp)
  %  ppp(i).xx
  %end
  %v0c = v0s/norm(v0s) + mp;
  %v1c = v1s/norm(v1s) + mp;
  %v2c = v2s/norm(v2s) + mp;
  v0c = v0s + mp;
  v1c = v1s + mp;
  v2c = v2s + mp;
  plot3(mp(1,1), mp(2,1), mp(3,1), 'b.', 'MarkerSize', 10); hold on;
  plot3([mp(1,1) v0c(1,1)], [mp(2,1) v0c(2,1)], [mp(3,1) v0c(3,1)], 'g-', 'LineWidth', 2); hold on;
  plot3([mp(1,1) v1c(1,1)], [mp(2,1) v1c(2,1)], [mp(3,1) v1c(3,1)],'g-', 'LineWidth', 2); hold on;
  plot3([mp(1,1) v2c(1,1)], [mp(2,1) v2c(2,1)], [mp(3,1) v2c(3,1)],'g-', 'LineWidth', 2); hold on;


  cc = [0;0;0];
  [point, dist, tt, dd] = distance(cc, ppp(1).xx, ppp(4).xx, ppp(3).xx, ppp(2).xx);
  plot3(point(1), point(2), point(3), 'k.', 'MarkerSize', 10);
  [point, dist, tt, dd] = distance(cc, ppp(5).xx, ppp(6).xx, ppp(7).xx, ppp(8).xx);
  plot3(point(1), point(2), point(3), 'k.', 'MarkerSize', 10);
  [point, dist, tt, dd] = distance(cc, ppp(1).xx, ppp(5).xx, ppp(8).xx, ppp(4).xx);
  plot3(point(1), point(2), point(3), 'k.', 'MarkerSize', 10);
  [point, dist, tt, dd] = distance(cc, ppp(2).xx, ppp(3).xx, ppp(7).xx, ppp(6).xx);
  plot3(point(1), point(2), point(3), 'k.', 'MarkerSize', 10);
  [point, dist, tt, dd] = distance(cc, ppp(1).xx, ppp(2).xx, ppp(6).xx, ppp(5).xx);
  plot3(point(1), point(2), point(3), 'k.', 'MarkerSize', 10);
  [point, dist, tt, dd] = distance(cc, ppp(3).xx, ppp(4).xx, ppp(8).xx, ppp(7).xx);
  plot3(point(1), point(2), point(3), 'k.', 'MarkerSize', 10);
  %axis('equal');
end

function [point, dist, tt, dd] = distance(x0, x1, x2, x3, x4)

  normal = cross(x2 - x1, x3 - x1);
  normal = normal/norm(normal);
  %normal

  dist = dot(normal, (x0 - x1));
  point = x0 - normal*dist;
  
  x1x0 = x1 - point;
  x2x0 = x2 - point;
  x2x1 = x2 - x1;
  tt(1) = -dot(x1x0, x2x1)/norm(x2x1)^2;
  dd(1) = norm(cross(-x1x0, -x2x0))/norm(x2x1);

  x1x0 = x2 - point;
  x2x0 = x3 - point;
  x2x1 = x3 - x2;
  tt(2) = -dot(x1x0, x2x1)/norm(x2x1)^2;
  dd(2) = norm(cross(-x1x0, -x2x0))/norm(x2x1);

  x1x0 = x3 - point;
  x2x0 = x4 - point;
  x2x1 = x4 - x3;
  tt(3) = -dot(x1x0, x2x1)/norm(x2x1)^2;
  dd(3) = norm(cross(-x1x0, -x2x0))/norm(x2x1);

  x1x0 = x4 - point;
  x2x0 = x1 - point;
  x2x1 = x1 - x4;
  tt(4) = -dot(x1x0, x2x1)/norm(x2x1)^2;
  dd(4) = norm(cross(-x1x0, -x2x0))/norm(x2x1);
end

function [verts] = createBoxVertices(mm, ee, v0, v1, v2)
  verts(1).xx = mm - ee(1)*v0 - ee(2)*v1 - ee(3)*v2;
  verts(2).xx = mm + ee(1)*v0 - ee(2)*v1 - ee(3)*v2;
  verts(3).xx = mm + ee(1)*v0 + ee(2)*v1 - ee(3)*v2;
  verts(4).xx = mm - ee(1)*v0 + ee(2)*v1 - ee(3)*v2;
  verts(5).xx = mm - ee(1)*v0 - ee(2)*v1 + ee(3)*v2;
  verts(6).xx = mm + ee(1)*v0 - ee(2)*v1 + ee(3)*v2;
  verts(7).xx = mm + ee(1)*v0 + ee(2)*v1 + ee(3)*v2;
  verts(8).xx = mm - ee(1)*v0 + ee(2)*v1 + ee(3)*v2;
end

function drawBoxFaces(verts, color)
  drawFace(verts(1), verts(4), verts(3), verts(2), color);
  drawFace(verts(5), verts(6), verts(7), verts(8), color);
  drawFace(verts(1), verts(5), verts(8), verts(4), color);
  drawFace(verts(2), verts(3), verts(7), verts(6), color);
  drawFace(verts(1), verts(2), verts(6), verts(5), color);
  drawFace(verts(3), verts(4), verts(8), verts(7), color);
end

function drawFace(v1, v2, v3, v4, color)
  h1 = patch([v1.xx(1,1) v2.xx(1,1) v3.xx(1,1)],...
             [v1.xx(2,1) v2.xx(2,1) v3.xx(2,1)],...
             [v1.xx(3,1) v2.xx(3,1) v3.xx(3,1)],...
             color, 'edgecolor', 'none');
  h2 = patch([v1.xx(1,1) v3.xx(1,1) v4.xx(1,1)],...
             [v1.xx(2,1) v3.xx(2,1) v4.xx(2,1)],...
             [v1.xx(3,1) v3.xx(3,1) v4.xx(3,1)],...
             color, 'edgecolor', 'none');
  plot3([v1.xx(1,1) v2.xx(1,1)],...
        [v1.xx(2,1) v2.xx(2,1)],...
        [v1.xx(3,1) v2.xx(3,1)],...
        'r');
  plot3([v2.xx(1,1) v3.xx(1,1)],...
        [v2.xx(2,1) v3.xx(2,1)],...
        [v2.xx(3,1) v3.xx(3,1)],...
        'g');
  plot3([v3.xx(1,1) v4.xx(1,1)],...
        [v3.xx(2,1) v4.xx(2,1)],...
        [v3.xx(3,1) v4.xx(3,1)],...
        'b');
  plot3([v4.xx(1,1) v1.xx(1,1)],...
        [v4.xx(2,1) v1.xx(2,1)],...
        [v4.xx(3,1) v1.xx(3,1)],...
        'm');
end
