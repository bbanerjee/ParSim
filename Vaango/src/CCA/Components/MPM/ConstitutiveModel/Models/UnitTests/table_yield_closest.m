function table_yield_closest()
  format long e
  K = 1.0e5;
  G = 1.0e5;
  sqrtKG = sqrt(1.5*K/G);
  zfac = -sqrt(3);
  rfac = sqrt(2)*sqrtKG;

  p = [-10   10  400  800  1600  3200 6400];
  q = [  0  100  500  600   700   800   900];

  z = p*zfac;
  rprime =  q*rfac;
  plot(z,rprime, 'b-', 'LineWidth', 2); hold on;
  axis equal

  point = [-2000, 4000]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'mx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'mo', 'LineWidth', 3, 'Markersize', 7);

  point = [2000, 4000]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'kx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'ko', 'LineWidth', 3, 'Markersize', 7);

  point = [-3000, 0]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'bx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'bo', 'LineWidth', 3, 'Markersize', 7);

  point = [-3000 1000]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'rx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'ro', 'LineWidth', 3, 'Markersize', 7);

  point = [3000 1000]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'gx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'go', 'LineWidth', 3, 'Markersize', 7);

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
