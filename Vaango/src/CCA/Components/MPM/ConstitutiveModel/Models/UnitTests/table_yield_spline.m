function table_yield_spline()
  %format long e
  K = 1.0e5;
  G = 1.0e5;
  sqrtKG = sqrt(1.5*K/G);
  zfac = -sqrt(3);
  rfac = sqrt(2)*sqrtKG;

  p = [-10   10  400  800  1600  3200 6400];
  q = [  0  100  500  600   700   800   900];

  z = p*zfac;
  rprime =  q*rfac;
  %plot(z,rprime, 'b-', 'LineWidth', 2); hold on;
  %axis equal

  %p = [1 2 3 4 5 6 7 8];
  %q = [-1 1 -1 1 -1 1 -1 1];
  %p = [1 2 3 4 5 6];
  %q = [-1 1 -1 1 -1 1];

  plot(p, q, 'b-', 'LineWidth', 2); hold on;
  %axis equal

  [xspline, yspline] = B_spline_quadratic(p, q);
  plot(xspline,yspline, 'r-', 'LineWidth', 2);

  [xspline, yspline] = B_spline_cubic(p, q);
  plot(xspline,yspline, 'g-', 'LineWidth', 2);

  [xspline, yspline] = spline_catmull_rom(p, q);
  plot(xspline,yspline, 'm-', 'LineWidth', 2);

  [xspline, yspline] = spline_cubic(p, q);
  plot(xspline,yspline, 'k-', 'LineWidth', 2);

  %point = [-2000, 4000]
  %point(1) = point(1)*zfac;
  %point(2) = point(2)*rfac;
  %plot(point(1), point(2), 'mx', 'LineWidth', 3, 'Markersize', 7);

  %[xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  %plot(xc, yc, 'mo', 'LineWidth', 3, 'Markersize', 7);

  %point = [2000, 4000]
  %point(1) = point(1)*zfac;
  %point(2) = point(2)*rfac;
  %plot(point(1), point(2), 'kx', 'LineWidth', 3, 'Markersize', 7);

  %[xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  %plot(xc, yc, 'ko', 'LineWidth', 3, 'Markersize', 7);

  %point = [-3000, 0]
  %point(1) = point(1)*zfac;
  %point(2) = point(2)*rfac;
  %plot(point(1), point(2), 'bx', 'LineWidth', 3, 'Markersize', 7);

  %[xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  %plot(xc, yc, 'bo', 'LineWidth', 3, 'Markersize', 7);

  %point = [-3000 1000]
  %point(1) = point(1)*zfac;
  %point(2) = point(2)*rfac;
  %plot(point(1), point(2), 'rx', 'LineWidth', 3, 'Markersize', 7);

  %[xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  %plot(xc, yc, 'ro', 'LineWidth', 3, 'Markersize', 7);

  %point = [3000 1000]
  %point(1) = point(1)*zfac;
  %point(2) = point(2)*rfac;
  %plot(point(1), point(2), 'gx', 'LineWidth', 3, 'Markersize', 7);

  %[xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  %plot(xc, yc, 'go', 'LineWidth', 3, 'Markersize', 7);

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

function [xspline, yspline] = B_spline_quadratic(xvals, yvals)

  n = length(xvals)-1;
  k = 2;

  tvals = linspace(0, 1, 10);

  count = 0;
  for j=1 : n - k + 1
    [length(xvals) j j+1 j+2]
    for i=1:length(tvals)
      t = tvals(i);
      Pj = [xvals(j) yvals(j)];
      Pj1 = [xvals(j+1) yvals(j+1)];
      Pj2 = [xvals(j+2) yvals(j+2)];
      P = quadraticUniformBSpline(n, j, t, Pj, Pj1, Pj2);
      count = count+1;
      xspline(count) = P(1);
      yspline(count) = P(2);
    end
  end

end

function [P] = quadraticUniformBSpline(n, j, t, Pk, Pk1, Pk2)

  A = [1 t t^2];
  if (j == 1)
   M = 0.5*[2 0 0; -4 4 0; 2 -3 1];
  elseif (j == n-1)
   M = 0.5*[1 1 0; -2 2 0; 1 -3 2];
  else
   M = 0.5*[1 1 0; -2 2 0; 1 -2 1];
  endif
  Pkx = [Pk(1); Pk1(1); Pk2(1)];
  Pky = [Pk(2); Pk1(2); Pk2(2)];
  Px = A*M*Pkx;
  Py = A*M*Pky;
  P = [Px Py];

end

function [xspline, yspline] = B_spline_cubic(xvals, yvals)

  n = length(xvals)-1;
  k = 3;

  tvals = linspace(0, 1, 10);

  count = 0;
  for j=1 : n - k + 1
    [length(xvals) j j+1 j+2 j+3]
    for i=1:length(tvals)
      t = tvals(i);
      Pj = [xvals(j) yvals(j)];
      Pj1 = [xvals(j+1) yvals(j+1)];
      Pj2 = [xvals(j+2) yvals(j+2)];
      Pj3 = [xvals(j+3) yvals(j+3)];
      P = cubicUniformBSpline(n, j, t, Pj, Pj1, Pj2, Pj3);
      count = count+1;
      xspline(count) = P(1);
      yspline(count) = P(2);
    end
  end

end

function [P] = cubicUniformBSpline(n, j, t, Pk, Pk1, Pk2, Pk3)

  mval = 4;
  A = [t^3 t^2 t 1];
  if (j == 1)
   M = [-1 7/4 -11/12 1/6; 3 -9/2 3/2 0; -3 3 0 0; 1 0 0 0];
  elseif (j == 2)
   M = [-1/4 7/12 -1/2 1/6; 3/4 -5/4 1/2 0; -3/4 1/4 1/2 0; 1/4 7/12 1/6 0];
  elseif (j == n-3)
   M = [-1/6 1/2 -7/12 1/4; 1/2 -1 1/2 0; -1/2 0 1/2 0; 1/6 2/3 1/6 0];
  elseif (j == n-2)
   M = [-1/6 11/12 -7/4 1; 1/2 -5/4 3/4 0; -1/2 -1/4 3/4 0; 1/6 7/12 1/4 0];
  else
   M = zeros(4);
   for ii=1:4
     ival = ii - 1;
     for jj=1:4
       jval = jj - 1;
       t1 = 1/factorial(mval-1);
       t2 = factorial(mval-1)/(factorial(ival)*factorial(mval - 1 - ival)) ;
       for kval= jval:mval-1
         t3 = (mval - (kval+1))^(ival);
         t4 = (-1)^(kval-jval);
         t5 = factorial(mval)/(factorial(kval-jval)*factorial(mval-kval+jval));
         %terms = [t1 t2 t3 t4 t5]
         M(ii,jj) = M(ii,jj) + t3*t4*t5;
       end
       M(ii,jj) = M(ii,jj)*t1*t2;
     end
   end
   M
  endif
  Pkx = [Pk(1); Pk1(1); Pk2(1); Pk3(1)];
  Pky = [Pk(2); Pk1(2); Pk2(2); Pk3(2)];
  Px = A*M*Pkx;
  Py = A*M*Pky;
  P = [Px Py];

end

function [xspline, yspline] = spline_catmull_rom(xx, yy)

  length(xx)
  xvals = [xx(1)-1.0e-6 xx xx(length(xx))+1.0e-6]
  yvals = [yy(1) yy yy(length(yy))]
 
  %xvals = xx;
  %yvals = yy;
  n = length(xvals)

  tvals = linspace(0, 1, 10);

  count = 0;
  tau = 0.2;
  for j=1:n-3
    [length(xvals) j-2 j-1 j j+1]
    for i=1:length(tvals)
      t = tvals(i);
      Pjm2 = [xvals(j) yvals(j)];
      Pjm1 = [xvals(j+1) yvals(j+1)];
      Pj = [xvals(j+2) yvals(j+2)];
      Pj1 = [xvals(j+3) yvals(j+3)];
      P = catmullRomSpline(t, tau, Pjm2, Pjm1, Pj, Pj1);
      count = count+1;
      xspline(count) = P(1);
      yspline(count) = P(2);
    end
  end

end

function [P] = catmullRomSpline(t, tau, Pkm2, Pkm1, Pk, Pk1)

  A = [1 t t^2 t^3];
  Pkx = [Pkm2(1); Pkm1(1); Pk(1); Pk1(1)];
  Pky = [Pkm2(2); Pkm1(2); Pk(2); Pk1(2)];
  M = [0 1 0 0;-tau 0 tau 0;2*tau tau-3 3-2*tau -tau;-tau 2-tau tau-2 tau];
  Px = A*M*Pkx;
  Py = A*M*Pky;
  P = [Px Py];

end

function [xspline, yspline] = spline_cubic(xx, yy)

  n = length(xx);
  tvals = linspace(0, 1, 10);
  [ax, ay] = cubicSplineCoefficients(xx, yy);
  count = 0;
  for j=1:n-1
    for i=1:length(tvals)
      t = tvals(i);
      count = count+1;
      xspline(count) = ax(j,1) + ax(j,2)*t + ax(j,3)*t^2 + ax(j,4)*t^3;
      yspline(count) = ay(j,1) + ay(j,2)*t + ay(j,3)*t^2 + ay(j,4)*t^3;
    end
  end 
end

function [ax, ay] = cubicSplineCoefficients(xx, yy)

  n = length(xx);

  [Dx, Dy] = cubicSplineDerivatives(xx, yy);

  a = xx(1:n-1);
  b = Dx(1:n-1)';
  c = 3*(xx(2:n) - xx(1:n-1)) - 2*Dx(1:n-1)' - Dx(2:n)';
  d = 2*(xx(1:n-1) - xx(2:n)) + Dx(1:n-1)' + Dx(2:n)';
  ax = [a' b' c' d'];
  [xx(2:n)' - (a' + b' + c' + d') Dx(2:n) - (b' + 2*c' + 3*d')]

  a = yy(1:n-1);
  b = Dy(1:n-1)';
  c = 3*(yy(2:n) - yy(1:n-1)) - 2*Dy(1:n-1)' - Dy(2:n)';
  d = 2*(yy(1:n-1) - yy(2:n)) + Dy(1:n-1)' + Dy(2:n)';
  ay = [a' b' c' d'];
  
end

function [Dx, Dy] = cubicSplineDerivatives(xx, yy)

  n = length(xx);
  lower = [ones(1,n-1) 0];
  upper = [0 ones(1,n-1)];
  diagonal = [2 ones(1,n-2)*4 2];
  A = spdiags([lower' diagonal' upper'],[-1 0 1], n, n);
  full(A)
  b = [(xx(2) - xx(1)) (xx(3:n) - xx(1:n-2)) (xx(n) - xx(n-1))]*3;
  b
  Dx = A \ b';
  Dx
  b = [(yy(2) - yy(1)) (yy(3:n) - yy(1:n-2)) (yy(n) - yy(n-1))]*3;
  Dy = A \ b';
end
