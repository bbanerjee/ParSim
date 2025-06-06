function table_yield_closest_cap()
 table_yield_closest_matlab();
 table_yield_closest_computed();
end

function table_yield_closest_computed()

  colors = [27,158,119;217,95,2;117,112,179;231,41,138;102,166,30;230,171,2;166,118,29;102,102,102]/255;
  
  idx = 1;
  ZRTable = [-17.3205 -173.205 0; 16.9741 -1.73205 0; 17.3205 0 0; 16.9741 1.73205 0; -17.3205 173.205 0; -692.82 866.025 0; -1385.64 1039.23 0; -2419.67 1168.48 0; -2510.7 1175.37 0; -2601.04 1173.06 0; -2689.99 1161.31 0; -2776.89 1139.65 0; -2861.07 1103.93 0; -2941.89 1059.23 0; -3018.73 1005.84 0; -3091.02 944.088 0; -3158.2 874.421 0; -3219.75 797.356 0; -3275.22 713.49 0; -3324.17 623.496 0; -3366.25 528.113 0; -3401.11 428.141 0; -3428.51 324.433 0; -3448.23 217.884 0; -3460.13 109.423 0; -3464.1 0 0; -3460.13 -109.423 0] ;
  plot(ZRTable(:,1), ZRTable(:,2), 'c.-', 'Color', colors(idx,:), 'LineWidth', 1); hold on;
  axis equal

  idx = idx+1;
  ZRPoint = [3464.1 6928.2 0];
  ZRSpline = [-355.07 519.615 0; -377.597 542.421 0; -400.142 564.649 0; -422.707 586.299 0; -445.291 607.372 0; -467.894 627.868 0; -490.517 647.787 0; -513.159 667.128 0; -535.82 685.892 0; -558.5 704.079 0; -581.199 721.688 0; -603.918 738.72 0; -626.656 755.174 0; -649.413 771.051 0; -672.19 786.351 0; -694.985 801.073 0; -717.8 815.219 0; -740.635 828.786 0; -763.488 841.777 0; -786.361 854.19 0; -809.253 866.025 0; -832.164 877.284 0; -855.094 887.965 0; -878.044 898.068 0; -901.013 907.595 0; -924.001 916.544 0; -947.008 924.915 0; -970.035 932.709 0; -993.081 939.926 0; -1016.15 946.566 0; -1039.23 952.628 0];
  ZRClose = [-664.895 781.451 0];
  %plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'b-', 'Color', colors(idx,:), 'LineWidth', 2); 
  ZRClose = [-638.5554741647053 794.8336959595669 0];
  %plot(ZRSpline(:,1), ZRSpline(:,2), 'r-', 'Color', colors(idx,:), 'LineWidth', 2);
  %idx = idx+1;
  plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'r-', 'Color', colors(idx,:), 'LineWidth', 2); 

  idx = idx+1;
  ZRPoint = [-3464.1 6928.2 0];
  ZRSpline = [-2645.515121710858 1167.183378107211 0; -2648.479136199274 1166.786177161588 0; -2651.440864821653 1166.377964042565 0; -2654.400307577993 1165.958738750141 0; -2657.357464468296 1165.528501284319 0; -2660.312335492561 1165.087251645096 0; -2663.264920650788 1164.634989832473 0; -2666.215219942977 1164.17171584645 0; -2669.163233369127 1163.697429687028 0; -2672.10896092924 1163.212131354205 0; -2675.052402623315 1162.715820847983 0; -2677.993558451352 1162.20849816836 0; -2680.932428413352 1161.690163315337 0; -2683.869012509313 1161.160816288915 0; -2686.803310739236 1160.620457089093 0; -2689.735323103121 1160.069085715871 0; -2692.665049600968 1159.506702169249 0; -2695.592490232778 1158.933306449227 0; -2698.517644998549 1158.348898555805 0; -2701.440513898283 1157.753478488983 0; -2704.361096931979 1157.147046248761 0; -2707.279394099636 1156.529601835139 0; -2710.195405401256 1155.901145248117 0; -2713.109130836838 1155.261676487696 0; -2716.020570406381 1154.611195553874 0; -2718.929724109888 1153.949702446653 0; -2721.836591947355 1153.277197166031 0; -2724.741173918786 1152.59367971201 0; -2727.643470024178 1151.899150084588 0; -2730.543480263532 1151.193608283767 0; -2733.441204636848 1150.477054309546 0];
  ZRClose = [-2651.440864821653 1166.377964042565 0];
  %plot(ZRSpline(:,1), ZRSpline(:,2), 'g-', 'Color', colors(idx,:), 'LineWidth', 2);
  %plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'g-', 'Color', colors(idx,:), 'LineWidth', 2); 

  idx = idx+1;
  ZRPoint = [5196.15 0 0];
  ZRSpline = [17.1473 -0.866025 0; 17.1585 -0.80829 0; 17.1689 -0.750555 0; 17.1785 -0.69282 0; 17.1873 -0.635085 0; 17.1954 -0.57735 0; 17.2027 -0.519615 0; 17.2093 -0.46188 0; 17.215 -0.404145 0; 17.22 -0.34641 0; 17.2243 -0.288675 0; 17.2277 -0.23094 0; 17.2304 -0.173205 0; 17.2324 -0.11547 0; 17.2335 -0.057735 0; 17.2339 0 0; 17.2335 0.057735 0; 17.2324 0.11547 0; 17.2304 0.173205 0; 17.2277 0.23094 0; 17.2243 0.288675 0; 17.22 0.34641 0; 17.215 0.404145 0; 17.2093 0.46188 0; 17.2027 0.519615 0; 17.1954 0.57735 0; 17.1873 0.635085 0; 17.1785 0.69282 0; 17.1689 0.750555 0; 17.1585 0.80829 0; 17.1473 0.866025 0];
  ZRClose = [17.23390553531033 0 0];
  %plot(ZRSpline(:,1), ZRSpline(:,2), 'b-', 'Color', colors(idx,:), 'LineWidth', 2);
  %plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'b-', 'Color', colors(idx,:), 'LineWidth', 2); 

  idx = idx+1;
  ZRPoint = [5196.15 1732.05 0];
  ZRSpline = [-0.173205 87.4686 0;-1.67258 93.474 0;-3.88441 100.059 0;-6.80869 107.223 0;-10.4454 114.966 0;-14.7946 123.288 0;-19.8562 132.19 0;-25.6303 141.671 0;-32.1168 151.731 0;-39.3158 162.371 0;-47.2273 173.59 0;-55.8511 185.388 0;-65.1875 197.766 0;-75.2362 210.722 0;-85.9975 224.258 0;-97.4712 238.373 0;-109.657 253.068 0;-122.556 268.342 0;-136.167 284.195 0;-150.49 300.627 0;-165.526 317.639 0;-181.275 335.23 0;-197.736 353.4 0;-214.909 372.149 0;-232.795 391.478 0;-251.393 411.386 0;-270.703 431.873 0;-290.726 452.94 0;-311.462 474.586 0;-332.91 496.811 0;-355.07 519.615 0];
  ZRClose = [-1.672583729842345 93.47397083224905 0];
  ZRClose = [-4.851258445554237 105.6207524480543 0];
  %plot(ZRSpline(:,1), ZRSpline(:,2), 'm-', 'Color', colors(idx,:), 'LineWidth', 2);
  plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'm-', 'Color', colors(idx,:), 'LineWidth', 2); 

  idx = idx+1;
  ZRPoint = [-5196.15 1732.05 0];
  ZRSpline = [-3247.485696309073 755.422716382653 0;-3249.3309787475 752.6237906946425 0;-3251.169026510297 749.8180554908673 0;-3252.999839597465 747.0055107713273 0;-3254.823418009004 744.1861565360225 0;-3256.639761744913 741.359992784953 0;-3258.448870805194 738.5270195181185 0;-3260.250745189846 735.6872367355195 0;-3262.045384898869 732.8406444371554 0;-3263.832789932262 729.9872426230268 0;-3265.612960290026 727.1270312931333 0;-3267.385895972162 724.2600104474751 0;-3269.151596978668 721.3861800860521 0;-3270.910063309545 718.5055402088644 0;-3272.661294964793 715.6180908159118 0;-3274.405291944412 712.7238319071946 0;-3276.142054248402 709.8227634827125 0;-3277.871581876762 706.9148855424656 0;-3279.593874829494 704.000198086454 0;-3281.308933106596 701.0787011146775 0;-3283.01675670807 698.1503946271364 0;-3284.717345633914 695.2152786238304 0;-3286.41069988413 692.2733531047596 0;-3288.096819458715 689.3246180699242 0;-3289.775704357672 686.3690735193237 0;-3291.447354581001 683.4067194529587 0;-3293.111770128699 680.4375558708289 0;-3294.768951000769 677.4615827729342 0;-3296.418897197209 674.4788001592749 0;-3298.061608718021 671.4892080298507 0;-3299.697085563203 668.4928063846618 0];
  ZRClose = [-3294.768951000769 677.4615827729342 0];
  %plot(ZRSpline(:,1), ZRSpline(:,2), 'k-', 'Color', colors(idx,:), 'LineWidth', 2);
  %plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'k-', 'Color', colors(idx,:), 'LineWidth', 2); 

  idx = idx+1;
  ZRPoint = [-5196.15 0 0];
  ZRSpline = [-3462.11 54.7115 0; -3462.24 51.0641 0; -3462.36 47.4166 0; -3462.47 43.7692 0; -3462.57 40.1218 0; -3462.67 36.4743 0; -3462.75 32.8269 0; -3462.83 29.1795 0; -3462.89 25.532 0; -3462.95 21.8846 0; -3463 18.2372 0; -3463.04 14.5897 0; -3463.07 10.9423 0; -3463.09 7.29487 0; -3463.1 3.64743 0; -3463.11 0 0; -3463.1 -3.64743 0; -3463.09 -7.29487 0; -3463.07 -10.9423 0; -3463.04 -14.5897 0; -3463 -18.2372 0; -3462.95 -21.8846 0; -3462.89 -25.532 0; -3462.83 -29.1795 0; -3462.75 -32.8269 0; -3462.67 -36.4743 0; -3462.57 -40.1218 0; -3462.47 -43.7692 0; -3462.36 -47.4166 0; -3462.24 -51.0641 0; -3462.11 -54.7115 0];
  ZRClose = [-3463.108025469086 0 0];
  %plot(ZRSpline(:,1), ZRSpline(:,2), 'k-', 'Color', colors(idx,:), 'LineWidth', 2);
  %plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'k-', 'Color', colors(idx,:), 'LineWidth', 2); 
end

function table_yield_closest_matlab()

  format long e
  K = 1.0e5;
  G = 1.0e5;
  sqrtKG = sqrt(1.5*K/G);
  zfac = -sqrt(3);
  rfac = sqrt(2)*sqrtKG;

  p = [-10   10  400  800  1600  3200 6400];
  q = [  0  100  500  600   700   800   900];

  R = 0.7;
  X = 3*2000;
  %X = 3*6400;
  %X = 3*10000;

  [p, q] = compute_cap(p, q, R, X);

  z = p*zfac;
  rprime =  q*rfac;
  nn = length(z)-1;
  for i=1:length(z)-1
    plot(z(i:i+1),rprime(i:i+1), '-.', 'Color', [i/nn 0 (1 - i/nn)], 'LineWidth', 3); hold on;
  end
  axis equal

  point = [-2000, 4000]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'rx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'rs', 'LineWidth', 3, 'Markersize', 7);
  [xc, yc, t] = closestPointNewtonSpline(point(1), point(2), z, rprime, [1 0 0])
  plot(xc, yc, 'rx', 'LineWidth', 3, 'Markersize', 7);

  point = [2000, 4000]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'gx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'gs', 'LineWidth', 3, 'Markersize', 7);
  [xc, yc, t] = closestPointNewtonSpline(point(1), point(2), z, rprime, [0 1 0])
  plot(xc, yc, 'gx', 'LineWidth', 3, 'Markersize', 7);

  point = [-3000, 0]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'bx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'bs', 'LineWidth', 3, 'Markersize', 7);
  [xc, yc, t] = closestPointNewtonSpline(point(1), point(2), z, rprime, [0 0 1])
  plot(xc, yc, 'bx', 'LineWidth', 3, 'Markersize', 7);

  point = [-3000 1000]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'mx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'ms', 'LineWidth', 3, 'Markersize', 7);
  [xc, yc, t] = closestPointNewtonSpline(point(1), point(2), z, rprime, [0.75 0 0.25])
  plot(xc, yc, 'mx', 'LineWidth', 3, 'Markersize', 7);


  point = [3000 1000]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'kx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'ks', 'LineWidth', 3, 'Markersize', 7);
  [xc, yc, t] = closestPointNewtonSpline(point(1), point(2), z, rprime, [0 0 0])
  plot(xc, yc, 'kx', 'LineWidth', 3, 'Markersize', 7);

  point = [3000 0]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'x', 'Color', [0.2 0.5 0.1], 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 's', 'Color', [0.2 0.5 0.1], 'LineWidth', 3, 'Markersize', 7);
  [xc, yc, t] = closestPointNewtonSpline(point(1), point(2), z, rprime, [0.2 0.5 0.1])
  plot(xc, yc, 'x', 'Color', [0.2 0.5 0.1], 'LineWidth', 3, 'Markersize', 7);

  point = [-2.65232499416893e+03 1.16625390832548e+03]
  plot(point(1), point(2), 'x', 'Color', [0.5 0.2 0.1], 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 's', 'Color', [0.5 0.2 0.1], 'LineWidth', 3, 'Markersize', 7);
  [xc, yc, t] = closestPointNewtonSpline(point(1), point(2), z, rprime, [0.2 0.5 0.1])
  plot(xc, yc, 'x', 'Color', [0.5 0.2 0.1], 'LineWidth', 3, 'Markersize', 7);

end

function [xc, yc, mindist] = closestPoint(xp, yp, xpoly, ypoly)

  xclose = [];
  yclose = [];
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

function [index] = closestSegment(xp, yp, xpoly, ypoly)

  index = 1;
  mindSq = 1.0e20;
  for i=1:length(xpoly)-1
   dSq = (xp - xpoly(i))^2 + (yp - ypoly(i))^2;
   if (dSq < mindSq) 
     mindSq = dSq;
     index = i;
   endif
  end
end

function [p_cap, q_cap] = compute_cap(p, q, R, X)
  p_min = min(p)
  p_max = X/3
  kappa = p_min + R*(p_max - p_min)
  
  count = 0;
  p_end = p(length(p));
  if (p_max > p(length(p)))
    p_start = p(length(p)-1);
    q_start = q(length(p)-1);
    q_end = q(length(p));
    dp = p_end - p_start;
    curr_p = p_end + dp;
    t = (curr_p - p_start)/dp;
    curr_q = (1 - t)*q_start + t*q_end;
    count = count+1;
    p_ext(count) = curr_p;
    q_ext(count) = curr_q;
    while (curr_p < p_max)
      curr_p = curr_p + dp;
      t = (curr_p - p_start)/dp;
      curr_q = (1 - t)*q_start + t*q_end;
      count = count+1;
      p_ext(count) = curr_p;
      q_ext(count) = curr_q;
    endwhile
    p = [p p_ext];
    q = [q q_ext];
  endif

  startp = 1;
  for i = 1:length(p)
    if (kappa > p(i)) 
      startp = i;
    endif 
  end
  pp = p(1:startp);
  qq = q(1:startp);

  count = 0;
  if (p_max > p(length(p)))
    p_start = p(length(p)-1);
    p_end = p(length(p));
    q_start = q(length(p)-1);
    q_end = q(length(p));
    dp = p_end - p_start;
    curr_p = p_start + dp;
    t = (curr_p - p_start)/dp;
    curr_q = (1 - t)*q_start + t*q_end;
    count = count+1;
    p_ext(count) = curr_p;
    q_ext(count) = curr_q;
    while (curr_p < p_max)
      curr_p = curr_p + dp;
      t = (curr_p - p_start)/dp;
      curr_q = (1 - t)*q_start + t*q_end;
      count = count+1;
      p_ext(count) = curr_p;
      q_ext(count) = curr_q;
    endwhile
    pp = [pp p_ext];
    qq = [qq q_ext];
  endif

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

function [B, T, N] = computeBSpline(loc, p1, p2, p3, t)

  spline_mat = [[1  1  0];[ -2  2  0];[ 1  -2  1]];
  if (loc == 'stt')
    spline_mat = [[2  0  0];[ -4  4  0];[ 2  -3  1]];
  elseif (loc == 'end')
    spline_mat = [[1  1  0];[ -2  2  0];[ 1  -3  2]];
  endif 

  A = [1 t t^2];
  dA = [0 1 2*t];
  ddA = [0 0 2];
  Px = [p1(1) p2(1) p3(1)];
  Py = [p1(2) p2(2) p3(2)];
  bx = dot(A, spline_mat * Px' * 0.5);
  by = dot(A, spline_mat * Py' * 0.5);
  tx = dot(dA, spline_mat * Px' * 0.5);
  ty = dot(dA, spline_mat * Py' * 0.5);
  nx = dot(ddA, spline_mat * Px' * 0.5);
  ny = dot(ddA, spline_mat * Py' * 0.5);
  B = [bx by];
  T = [tx ty];
  N = [nx ny];
end

function [xc, yc, t] = closestPointNewtonSpline(px, py, z, rprime, color)

  [index] = closestSegment(px, py, z, rprime);
  %plot([z(index) z(index+1)], [rprime(index) rprime(index+1)], '-', 'Color', color, 'LineWidth', 4, 'Markersize', 7);
  if (index < 3)
    loc = 'stt';
    p1 = [z(1) rprime(1)];
    p2 = [z(2) rprime(2)];
    p3 = [z(3) rprime(3)];
  elseif (index > length(z)-2)
    loc = 'end';
    p1 = [z(length(z)-2) rprime(length(z)-2)];
    p2 = [z(length(z)-1) rprime(length(z)-1)];
    p3 = [z(length(z)) rprime(length(z))];
  else
    loc = 'mid';
    p1 = [z(index-1) rprime(index-1)];
    p2 = [z(index) rprime(index)];
    p3 = [z(index+1) rprime(index+1)];
  endif
  %plot([p1(1) p2(1) p3(1)],[p1(2) p2(2) p3(2)], 'g-');
  t = linspace(0, 1, 20);
  for i=1:length(t)
    [B, T, N] = computeBSpline(loc, p1, p2, p3, t(i));
    T_hat = T/norm(T);
    N_hat = N/norm(N);
    plot(B(1), B(2), 'r.');
    %plot([B(1) B(1) + T_hat(1)*10], [B(2) B(2) + T_hat(2)*10], 'b-');
    %plot([B(1) B(1) + N_hat(1)*10], [B(2) B(2) + N_hat(2)*10], 'g-');
  end
  [xc, yc, t] = closestPointNewton(px, py, loc, p1, p2, p3);
end

function [xc, yc, t] = closestPointNewton(xp, yp, loc, p1, p2, p3)

  p = [xp yp];
  t = 0.5;
  f = 1;
  nn = 0;
  while (abs(f) > 1.0e-10) 
    nn = nn+1;
    [B, T, N] = computeBSpline(loc, p1, p2, p3, t);
    f = dot(p - B, T);
    fp = -dot(T, T) + dot(p - B, N);
    t = t - f/fp;
    %[nn f fp t]
    if (nn > 20) 
      [f t fp]
      break;
    end
  endwhile

  [B, T, N] = computeBSpline(loc, p1, p2, p3, t);
  xc = B(1);
  yc = B(2);
  
end
