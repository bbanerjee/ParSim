function table_yield()
  p = [-10   10  200  300  400  800  1600  3200 6400];
  q = [  0  100  200  300  500  600   700   800   900];
  plot(p,q); hold on;

  [p, q] = hull(p, q);
  plot(p, q,'g-', 'LineWidth', 2);
  [p' q']

  p_closed = [p p(length(p)) p(1)];
  q_closed = [q q(1) q(1)];
  plot(p_closed, q_closed,'m-');

  [p_cen, q_cen] = calcCentroid(p_closed, q_closed);
  plot(p_cen, q_cen, 'x');

  p_scale = p_cen + 1.5*(p_closed - p_cen);
  q_scale = q_cen + 1.5*(q_closed - q_cen);
  plot(p_scale, q_scale,'g-');

  p_trans = p_scale;
  q_trans = q_scale - min(q_scale);
  plot(p_trans, q_trans, 'r-');

  [xnormal, ynormal, p_ext, q_ext] = calcNormals(p, q);
  plot(p_ext, q_ext, 'r-', 'Linewidth', 3);

  p_normal = p + xnormal*900;
  q_normal = q + ynormal*900;
  for i=1:length(p)
    plot([p(i) p_normal(i)],[q(i) q_normal(i)], 'b-', 'Linewidth', 1);
  end
  
  plot(p_normal, q_normal,'b-', 'LineWidth', 2);
  axis equal
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
  [xx' yy']

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
  [xgrad' ygrad']

  gradx = gradient(xx(:));
  grady = gradient(yy(:));
  dr = [gradx grady];
  size(dr);
  T = [];
  for i=1:size(dr)(1)
   T(i,:) = dr(i,:)/(sqrt(dr(i,1)^2 + dr(i,2)^2));
  end
  T

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
  N
 
  xnormal = -xnormal;
  ynormal = -ynormal;
  [xnormal' ynormal']
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
