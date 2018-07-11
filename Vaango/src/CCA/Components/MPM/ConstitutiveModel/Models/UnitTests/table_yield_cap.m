function table_yield_cap()

  close all;

  fig = figure()

  p_orig = [-10   10  200  300  400  800  1600  3200 6400 6720 6752];
  q_orig = [  0  100  200  300  500  600   700   800   900 910 911];
  plot(p_orig, q_orig, 'k.', 'Markersize', 10); hold on;

  [p, q] = hull(p_orig, q_orig);
  plot(p, q, '-', 'Color', [216,179,101]/255, 'LineWidth', 2);

  R = 0.7;
  X = 3*2000;
  %X = 3*6400;
  %X = 3*10000;

  [p, q] = compute_cap(p, q, R, X);
  %p_q = [p' q']
  %plot(p, q, '-', 'Color', [116,179,101]/255, 'LineWidth', 2);
  
  %axis equal
  %p_q = [p' q']
  %print -dpdf tabular_yield_hull.pdf

  %plot_stress_states(p_orig, q_orig, R, X, fig)

  q_max = [max(q_orig) max(q)]

  plot_table_yield(p, q);
end

function [p_cap, q_cap] = compute_cap(p, q, R, X)
  p_min = min(p)
  p_max = X/3
  kappa = p_min + R*(p_max - p_min)
  
  startp = 1;
  for i = 1:length(p)
    if (kappa > p(i)) 
      startp = i;
    endif 
  end
  pp = p(1:startp);
  qq = q(1:startp);

  theta = linspace(-5*pi/180, pi/2, 20);
  theta = fliplr(theta);

  a = p_max - kappa;

  p_cap = kappa + a*cos(theta);
  for i=1:length(theta)
    b(i) = computeEllipseHeight(p, q, p_cap(i));
  end
  q_cap = b.*sin(theta);
  bb = b;

  p_cap = [pp p_cap];
  q_cap = [qq q_cap];

end

function [q_cap] = computeEllipseHeight(p, q, p_cap)
  startp = 1;
  endp = 2;
  for i = 1:length(p)
    if (p_cap > p(i)) 
      startp = i;
      endp = i+1;
    endif 
  end
  if (startp == length(p))
    startp = startp - 1;
    endp = endp - 1;
  endif
  t = (p_cap - p(startp))/(p(endp) - p(startp));
  q_cap = (1 - t)*q(startp) + t*q(endp);
  %tt = [startp endp p(startp) p(endp) t]
end 

function plot_stress_states(plist, qlist, R, X, fig)

  p_min = min(plist)
  p_max = X/3;
  kappa = p_min + R*(p_max - p_min)
  a = p_max - kappa;

  p = -300; q = 1000;
  plot(p, q, 'ro');

  p = -2; q = 39;
  plot(p, q, 'go');

  p = 1000; q = 605;
  plot(p, q, 'bo');

  p = 1000; q = 635;
  plot(p, q, 'co');

  p = 7000; q = 1000;
  plot(p, q, 'mo');

  p = 2000; q = 0;
  plot(p, q, 'yo');

  %p = 1700; % for p_max = 2000
  p = 9700; % for p_max = 2000
  [q_val] = computeEllipseHeight(plist, qlist, p)
  theta = acos((p - kappa)/a);
  q = q_val*sin(theta);
  [p, q]
  plot(p, q, 'ko');
end

function plot_table_yield(p, q)

  grid minor

  p_closed = [p p(1)];
  q_closed = [q q(1)];
  plot(p_closed, q_closed,'g--');

  [p_cen, q_cen] = calcCentroid(p_closed, q_closed);
  plot(p_cen, q_cen, 'x');

  [xnormal, ynormal, p_ext, q_ext] = calcNormals(p, q);
  %plot(p_ext, q_ext, 'r-', 'Linewidth', 3);

  p_normal = p + xnormal*900;
  q_normal = q + ynormal*900;
  for i=1:length(p)
    plot([p(i) p_normal(i)],[q(i) q_normal(i)], 'b-', 'Linewidth', 1);
  end
  
  plot(p_normal, q_normal,'b-', 'LineWidth', 2);
  axis equal

  point = [0, 0];
  plot(point(1), point(2), 'x', 'Color', [116,179,101]/255, 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), p, q);
  plot(xc, yc, 'o', 'Color', [116,179,101]/255, 'LineWidth', 3, 'Markersize', 7);

  point = [-2000, 0];
  plot(point(1), point(2), 'gx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), p, q);
  plot(xc, yc, 'go', 'LineWidth', 3, 'Markersize', 7);

  point = [-300, 1000];
  plot(point(1), point(2), 'mx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), p, q);
  plot(xc, yc, 'mo', 'LineWidth', 3, 'Markersize', 7);

  %p_new = p + xnormal*r;
  %q_new = q + ynormal*r;
  %plot(p_new, q_new,'m-', 'LineWidth', 2);

  point = [-2000, 4000];
  plot(point(1), point(2), 'cx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), p, q);
  plot(xc, yc, 'co', 'LineWidth', 3, 'Markersize', 7);

  %p_new = p + xnormal*r;
  %q_new = q + ynormal*r;
  %plot(p_new, q_new,'m-', 'LineWidth', 2);

  point = [2000, 4000];
  plot(point(1), point(2), 'kx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), p, q);
  plot(xc, yc, 'ko', 'LineWidth', 3, 'Markersize', 7);

  %p_new = p + xnormal*r;
  %q_new = q + ynormal*r;
  %plot(p_new, q_new,'k-', 'LineWidth', 2);

  point = [3000, 0];
  plot(point(1), point(2), 'bx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), p, q);
  plot(xc, yc, 'bo', 'LineWidth', 3, 'Markersize', 7);

  %p_new = p + xnormal*r;
  %q_new = q + ynormal*r;
  %plot(p_new, q_new,'b-', 'LineWidth', 2);

  point = [3000 1000];
  plot(point(1), point(2), 'rx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), p, q);
  plot(xc, yc, 'ro', 'LineWidth', 3, 'Markersize', 7);
  closest = [xc, yc, r]

  %p_new = p + xnormal*r;
  %q_new = q + ynormal*r;
  %plot(p_new, q_new,'r-', 'LineWidth', 2);
end

function [x_cen, y_cen] = calcCentroid(x, y) 
  area = 0;
  x_cen = 0;
  y_cen = 0;
  for i=1:length(x)-1
    area_term = x(i)*y(i+1) - x(i+1)*y(i);
    area = area + area_term;
    x_cen = x_cen + (x(i) + x(i+1))*area_term;
    y_cen = y_cen + (y(i) + y(i+1))*area_term;
  end
  area = area/2.0;
  x_cen = x_cen/(6*area);
  y_cen = y_cen/(6*area);
end

function [xnormal, ynormal, xx, yy] = calcNormals(x, y)

  xx0 = x(3);
  yy0 = -y(3);
  xx1 = x(2);
  yy1 = -y(2);
  n = length(x);
  t = 1.1;
  xxn1 = (1 - t)*x(n-1) + t*x(n);
  yyn1 = (1 - t)*y(n-1) + t*y(n);
  xxn2 = (1 - t)*x(n) + t*xxn1;
  yyn2 = (1 - t)*y(n) + t*yyn1;
  xx = [xx0 xx1 x xxn1 xxn2];
  yy = [yy0 yy1 y yyn1 yyn2];
  %[xx' yy']

  xgrad = [];
  ygrad = [];
  for i=2:length(xx)-1
    tminus = 0.99;
    tplus = 0.01;
    xminus = (1 - tminus)*xx(i-1) + tminus*xx(i);
    xplus  = (1 - tplus)*xx(i) + tplus*xx(i+1);
    yminus = (1 - tminus)*yy(i-1) + tminus*yy(i);
    yplus  = (1 - tplus)*yy(i) + tplus*yy(i+1);
    xgrad(i-1) = (xplus - xminus)/0.02;
    ygrad(i-1) = (yplus - yminus)/0.02;
  end
  %[xgrad' ygrad']

  for i=1:length(xgrad)
    len = sqrt(xgrad(i)^2 + ygrad(i)^2);
    xgrad(i) = xgrad(i)/len;
    ygrad(i) = ygrad(i)/len;
  end
  %[xgrad' ygrad']

  gradx = gradient(xx(:));
  grady = gradient(yy(:));
  dr = [gradx grady];
  size(dr);
  T = [];
  for i=1:size(dr)(1)
   T(i,:) = dr(i,:)/(sqrt(dr(i,1)^2 + dr(i,2)^2));
  end
  %T

  %xnormal = xgrad(2:length(xgrad)-1)
  %ynormal = ygrad(2:length(ygrad)-1)

  xnormal = [];
  ynormal = [];
  for i=2:length(xgrad)-1
    tminus = 0.99;
    tplus = 0.01;
    xminus = (1 - tminus)*xgrad(i-1) + tminus*xgrad(i);
    xplus  = (1 - tplus)*xgrad(i) + tplus*xgrad(i+1);
    yminus = (1 - tminus)*ygrad(i-1) + tminus*ygrad(i);
    yplus  = (1 - tplus)*ygrad(i) + tplus*ygrad(i+1);
    xnormal(i-1) = (xplus - xminus)/0.02;
    ynormal(i-1) = (yplus - yminus)/0.02;
  end
  %[xnormal' ynormal']

  for i=1:length(xnormal)
    len = sqrt(xnormal(i)^2 + ynormal(i)^2);
    xnormal(i) = xnormal(i)/len;
    ynormal(i) = ynormal(i)/len;
  end

  ggradx = gradient(gradx);
  ggrady = gradient(grady);

  Tgradx = gradient(T(:,1));
  Tgrady = gradient(T(:,2));
  ddr = [Tgradx Tgrady];

  %size(ddr);
  N = [];
  for i=1:size(ddr)(1)
   N(i,:) = ddr(i,:)/(sqrt(ddr(i,1)^2 + ddr(i,2)^2));
  end
  %N
 
  xnormal = -xnormal;
  ynormal = -ynormal;
  normals = [xnormal' ynormal'];
end

function [val] = ccw(p1x, p1y, p2x, p2y, p3x, p3y)
  val = (p2x - p1x)*(p3y - p1y) - (p2y - p1y)*(p3x - p1x);
end

function [hullx, hully] = hull(x, y)
  n = length(x);
  
  lowerx = [];
  lowery = [];
  for i=1:n
    p3x = x(i);
    p3y = y(i);
    k = length(lowerx);
    lower_val = -1;
    while (k > 1 && lower_val <= 0) 
      p1x = lowerx(k-1);
      p1y = lowery(k-1);
      p2x = lowerx(k);
      p2y = lowery(k);
      points = [k-1 k i];
      lower_val = ccw(p1x, p1y, p2x, p2y, p3x, p3y);
      if (lower_val <= 0) 
        lowerx = lowerx(1:k-1);
        lowery = lowery(1:k-1);
        k = length(lowerx);
      endif
    endwhile
    lowerx(k+1) = x(i);
    lowery(k+1) = y(i);
    [lowerx' lowery'];
    %plot(lowerx, lowery, 'b--');
  end

  upperx = [];
  uppery = [];
  for i=n:-1:1
    p3x = x(i);
    p3y = y(i);
    k = length(upperx);
    upper_val = -1;
    while (k > 1 && upper_val <=0 )
      p1x = upperx(k-1);
      p1y = uppery(k-1);
      p2x = upperx(k);
      p2y = uppery(k);
      points = [k-1 k i];
      upper_val = ccw(p1x, p1y, p2x, p2y, p3x, p3y);
      if (upper_val <= 0) 
        upperx = upperx(1:k-1);
        uppery = uppery(1:k-1);
        k = length(upperx);
      endif
    endwhile
    upperx(k+1) = x(i);
    uppery(k+1) = y(i);
    [upperx' uppery'];
    %plot(upperx, uppery, 'b--');
  end

  lowerx = lowerx(1:length(lowerx)-1);
  upperx = upperx(1:length(upperx)-1);
  lowery = lowery(1:length(lowery)-1);
  uppery = uppery(1:length(uppery)-1);

  hullx = [lowerx upperx];
  hully = [lowery uppery];

  combined = [hullx' hully'];
  combined = sortrows(combined);
  hullx = combined(:,1)';
  hully = combined(:,2)';
  
end

function [xc, yc, mindist] = closestPoint(xp, yp, xpoly, ypoly)

  xclose = [];
  xclose = [];
  npts = 0;
  for i = 2:length(xpoly)
    xstart = xpoly(i-1);
    ystart = ypoly(i-1);
    xnext = xpoly(i);
    ynext = ypoly(i);
    xm = xnext - xstart;
    ym = ynext - ystart;
    xn = xp - xstart;
    yn = yp - ystart;
    if ((xm^2 + ym^2) < 1.0e-12)
      npts = npts + 1;
      xclose(npts) = xstart;   
      yclose(npts) = ystart;   
    else
      t0 = (xm*xn + ym*yn)/(xm*xm + ym*ym);
      if (t0 <= 0.0)
        npts = npts + 1;
        xclose(npts) = xstart;   
        yclose(npts) = ystart;   
      elseif (t0 >= 1.0)
        npts = npts + 1;
        xclose(npts) = xnext;   
        yclose(npts) = ynext;   
      else
        xmin = xm * t0 + xstart;
        ymin = ym * t0 + ystart;
        npts = npts + 1;
        xclose(npts) = xmin;   
        yclose(npts) = ymin;   
      endif
    endif
  end
  mindSq = 1.0e20;
  for i=1:length(xclose)
   dSq = (xp - xclose(i))^2 + (yp - yclose(i))^2;
   if (dSq < mindSq) 
     mindSq = dSq;
     xc = xclose(i);
     yc = yclose(i);
   endif
  end
  mindist = sqrt(mindSq);
end
