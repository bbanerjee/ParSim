function table_yield_closest()
 %table_yield_closest_matlab();
 table_yield_closest_computed();
end

function table_yield_closest_computed()

  idx = 1;
  colors = [27,158,119;217,95,2;117,112,179;231,41,138;102,166,30;230,171,2;166,118,29;102,102,102]/255;
  ZRTable = [-692.82 -866.025 0; -17.3205 -173.205 0; 17.3205 0 0; -17.3205 173.205 0; -692.82 866.025 0; -1385.64 1039.23 0; -2771.28 1212.44 0; -5542.56 1385.64 0; -11085.1 1558.85 0; -11639.4 1576.17 0; -11694.8 1577.9 0] ;
  plot(ZRTable(:,1), ZRTable(:,2), 'c.-', 'Color', colors(idx,:), 'LineWidth', 1); hold on;
  axis equal

  idx = idx+1;
  ZRPoint = [3464.1 6928.2 0];
  ZRSpline = [-355.07 519.615 0; -377.597 542.421 0; -400.142 564.649 0; -422.707 586.299 0; -445.291 607.372 0; -467.894 627.868 0; -490.517 647.787 0; -513.159 667.128 0; -535.82 685.892 0; -558.5 704.079 0; -581.199 721.688 0; -603.918 738.72 0; -626.656 755.174 0; -649.413 771.051 0; -672.19 786.351 0; -694.985 801.073 0; -717.8 815.219 0; -740.635 828.786 0; -763.488 841.777 0; -786.361 854.19 0; -809.253 866.025 0; -832.164 877.284 0; -855.094 887.965 0; -878.044 898.068 0; -901.013 907.595 0; -924.001 916.544 0; -947.008 924.915 0; -970.035 932.709 0; -993.081 939.926 0; -1016.15 946.566 0] ;
  ZRClose = [-664.895 781.451 0];
  plot(ZRSpline(:,1), ZRSpline(:,2), 'r-', 'Color', colors(idx,:), 'LineWidth', 2);
  plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'r-', 'Color', colors(idx,:), 'LineWidth', 2); 

  idx = idx+1;
  ZRPoint = [-3464.1 6928.2 0];
  ZRSpline = [-2078.46 1125.83 0; -2125.42 1131.61 0; -2173.92 1137.38 0; -2223.95 1143.15 0; -2275.53 1148.93 0; -2328.65 1154.7 0; -2383.3 1160.47 0; -2439.5 1166.25 0; -2497.23 1172.02 0; -2556.51 1177.79 0; -2617.32 1183.57 0; -2679.68 1189.34 0; -2743.57 1195.12 0; -2809 1200.89 0; -2875.97 1206.66 0; -2944.49 1212.44 0; -3014.54 1218.21 0; -3086.13 1223.98 0; -3159.26 1229.76 0; -3233.93 1235.53 0; -3310.14 1241.3 0; -3387.89 1247.08 0; -3467.18 1252.85 0; -3548.01 1258.62 0; -3630.38 1264.4 0; -3714.29 1270.17 0; -3799.73 1275.94 0; -3886.72 1281.72 0; -3975.25 1287.49 0; -4065.32 1293.26 0] ;
  ZRClose = [-3839.08 1278.56 0];
  plot(ZRSpline(:,1), ZRSpline(:,2), 'g-', 'Color', colors(idx,:), 'LineWidth', 2);
  plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'g-', 'Color', colors(idx,:), 'LineWidth', 2); 

  idx = idx+1;
  ZRPoint = [5196.15 0 0];
  ZRSpline = [0 -86.6025 0; 1.11621 -80.829 0; 2.15544 -75.0555 0; 3.11769 -69.282 0; 4.00296 -63.5085 0; 4.81125 -57.735 0; 5.54256 -51.9615 0; 6.19689 -46.188 0; 6.77424 -40.4145 0; 7.27461 -34.641 0; 7.698 -28.8675 0; 8.04441 -23.094 0; 8.31384 -17.3205 0; 8.50629 -11.547 0; 8.62176 -5.7735 0; 8.66025 0 0; 8.62176 5.7735 0; 8.50629 11.547 0; 8.31384 17.3205 0; 8.04441 23.094 0; 7.698 28.8675 0; 7.27461 34.641 0; 6.77424 40.4145 0; 6.19689 46.188 0; 5.54256 51.9615 0; 4.81125 57.735 0; 4.00296 63.5085 0; 3.11769 69.282 0; 2.15544 75.0555 0; 1.11621 80.829 0] ;
  ZRClose = [8.66025 0 0];
  plot(ZRSpline(:,1), ZRSpline(:,2), 'b-', 'Color', colors(idx,:), 'LineWidth', 2);
  plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'b-', 'Color', colors(idx,:), 'LineWidth', 2); 

  idx = idx+1;
  ZRPoint = [5196.15 1732.05 0];
  ZRSpline = [0 86.6025 0; -1.51073 92.6647 0; -3.73353 99.3042 0; -6.6684 106.521 0; -10.3153 114.315 0; -14.6743 122.687 0; -19.7454 131.636 0; -25.5285 141.162 0; -32.0237 151.266 0; -39.231 161.947 0; -47.1503 173.205 0; -55.7817 185.041 0; -65.1251 197.454 0; -75.1806 210.444 0; -85.9482 224.012 0; -97.4279 238.157 0; -109.62 252.879 0; -122.523 268.179 0; -136.139 284.056 0; -150.467 300.511 0; -165.507 317.543 0; -181.259 335.152 0; -197.723 353.338 0; -214.899 372.102 0; -232.788 391.443 0; -251.388 411.362 0; -270.7 431.858 0; -290.725 452.931 0; -311.461 474.582 0; -332.91 496.81 0] ;
  ZRClose = [-1.51073 92.6647 0];
  plot(ZRSpline(:,1), ZRSpline(:,2), 'm-', 'Color', colors(idx,:), 'LineWidth', 2);
  plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'm-', 'Color', colors(idx,:), 'LineWidth', 2); 

  idx = idx+1;
  ZRPoint = [-5196.15 1732.05 0];
  ZRSpline = [-4156.92 1299.04 0; -4250.84 1304.81 0; -4347.83 1310.59 0; -4447.91 1316.36 0; -4551.06 1322.13 0; -4657.29 1327.91 0; -4766.6 1333.68 0; -4878.99 1339.45 0; -4994.46 1345.23 0; -5113.01 1351 0; -5234.64 1356.77 0; -5359.35 1362.55 0; -5487.14 1368.32 0; -5618 1374.09 0; -5751.95 1379.87 0; -5888.97 1385.64 0; -6029.08 1391.41 0; -6172.26 1397.19 0; -6318.52 1402.96 0; -6467.86 1408.73 0; -6620.28 1414.51 0; -6775.78 1420.28 0; -6934.36 1426.06 0; -7096.02 1431.83 0; -7260.76 1437.6 0; -7428.57 1443.38 0; -7599.47 1449.15 0; -7773.44 1454.92 0; -7950.5 1460.7 0; -8130.63 1466.47 0] ;
  ZRClose = [-5214.01 1355.79 0];
  plot(ZRSpline(:,1), ZRSpline(:,2), 'k-', 'Color', colors(idx,:), 'LineWidth', 2);
  plot([ZRPoint(1) ZRClose(1)], [ZRPoint(2) ZRClose(2)], 'k-', 'Color', colors(idx,:), 'LineWidth', 2); 
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

  z = p*zfac;
  rprime =  q*rfac;
  nn = length(z)-1;
  for i=1:length(z)-1
    plot(z(i:i+1),rprime(i:i+1), '-s', 'Color', [i/nn 0 (1 - i/nn)], 'LineWidth', 3); hold on;
  end
  axis equal

  point = [-2000, 4000]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'rx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'ro', 'LineWidth', 3, 'Markersize', 7);

  point = [2000, 4000]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'gx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'go', 'LineWidth', 3, 'Markersize', 7);

  point = [-3000, 0]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'bx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'bo', 'LineWidth', 3, 'Markersize', 7);

  point = [-3000 1000]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'mx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'mo', 'LineWidth', 3, 'Markersize', 7);

  point = [3000 1000]
  point(1) = point(1)*zfac;
  point(2) = point(2)*rfac;
  plot(point(1), point(2), 'kx', 'LineWidth', 3, 'Markersize', 7);

  [xc, yc, r] = closestPoint(point(1), point(2), z, rprime)
  plot(xc, yc, 'ko', 'LineWidth', 3, 'Markersize', 7);

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
