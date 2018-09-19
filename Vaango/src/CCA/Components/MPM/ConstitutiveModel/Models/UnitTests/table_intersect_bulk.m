function table_intersect_bulk()

  testdata()
  %drysand()

end

function testdata()

  color_00 = [18 150 155] ./ 255;
  color_01 = [94 250 81] ./ 255;
  color_02 = [12 195 82] ./ 255;

  poly_00 = [[0, 0];[0.1, 1000]];
  poly_01 = [[0.5 0];[0.55 1100];[0.75 1500];[0.85 2000]];
  poly_02 = [[1.0 0];[1.45 2000];[1.75 3000]];

  fig1 = figure();
  plot(poly_00(:,1), poly_00(:,2), '-', 'color', color_00, 'linewidth', 2); hold on;
  plot(poly_01(:,1), poly_01(:,2), '-', 'color', color_01, 'linewidth', 2); hold on;
  plot(poly_02(:,1), poly_02(:,2), '-', 'color', color_02, 'linewidth', 2); hold on;
  extendPolyLoading(poly_00, color_00, 1.75);
  extendPolyLoading(poly_01, color_01, 1.75);
  extendPolyLoading(poly_02, color_02, 1.75);

  p = 0.0;
  dp = 1.0e-6;
  seg = [[-1 p - dp];[1.75 p - dp]];
  plot(seg(:,1), seg(:,2), 'r--', 'linewidth', 2);
  [p_00_lo, t_00_lo] = lineIntersectSegLinear(seg, poly_00)
  [p_01_lo, t_01_lo] = lineIntersectSegLinear(seg, poly_01)
  [p_02_lo, t_02_lo] = lineIntersectSegLinear(seg, poly_02);

  seg = [[-1 p + dp];[1.75 p + dp]];
  [p_00_hi, t_00_hi] = lineIntersectSegLinear(seg, poly_00)
  [p_01_hi, t_01_hi] = lineIntersectSegLinear(seg, poly_01)
  [p_02_hi, t_02_hi] = lineIntersectSegLinear(seg, poly_02);

  fig2 = figure();
  eps_p = [0, 0.5, 1.0];
  eps_v_lo = [p_00_lo(1), p_01_lo(1), p_02_lo(1)];
  eps_v_hi = [p_00_hi(1), p_01_hi(1), p_02_hi(1)];
  plot(eps_v_lo, eps_p, '-x', 'color', [1 0 0], 'linewidth', 2); hold on;
  plot(eps_v_hi, eps_p, '-x', 'color', [0 0 1], 'linewidth', 2); hold on;

  seg = [[0.06 -1];[0.06 1]];
  poly_lo = [eps_v_lo' eps_p'];
  poly_hi = [eps_v_hi' eps_p'];
  [eps_lo, t_lo] = lineIntersectSegLinear(seg, poly_lo)
  [eps_hi, t_hi] = lineIntersectSegLinear(seg, poly_hi)
  plot(eps_lo(1), eps_lo(2), 'kx', 'linewidth', 2, 'markersize', 8);
  plot(eps_hi(1), eps_hi(2), 'kx', 'linewidth', 2, 'markersize', 8);

  dp = 2 * dp;
  deps_v = eps_hi(2) - eps_lo(2);
  K = dp/deps_v
end

function drysand()

  color_00 = [18 150 155] ./ 255;
  color_09 = [94 250 81] ./ 255;
  color_18 = [12 195 82] ./ 255;
  color_27 = [8 180 238] ./ 255;
  color_36 = [1 17 181] ./ 255;
  color_45 = [251 111 66] ./ 255;

  poly_00 = load('Sand_load_00.csv');
  poly_09 = load('Sand_unload_09.csv');
  poly_18 = load('Sand_unload_18.csv');
  poly_27 = load('Sand_unload_27.csv');
  poly_36 = load('Sand_unload_36.csv');
  poly_45 = load('Sand_unload_45.csv');
  [poly_00] = changeUnits(poly_00);
  [poly_09] = changeUnits(poly_09);
  [poly_18] = changeUnits(poly_18);
  [poly_27] = changeUnits(poly_27);
  [poly_36] = changeUnits(poly_36);
  [poly_45] = changeUnits(poly_45);

  fig1 = figure();
  plot(poly_00(:,1), poly_00(:,2), '-x', 'color', color_00, 'linewidth', 2); hold on;
  plot(poly_09(:,1), poly_09(:,2), '-', 'color', color_09, 'linewidth', 2); hold on;
  plot(poly_18(:,1), poly_18(:,2), '-', 'color', color_18, 'linewidth', 2); hold on;
  plot(poly_27(:,1), poly_27(:,2), '-', 'color', color_27, 'linewidth', 2); hold on;
  plot(poly_36(:,1), poly_36(:,2), '-', 'color', color_36, 'linewidth', 2); hold on;
  plot(poly_45(:,1), poly_45(:,2), '-', 'color', color_45, 'linewidth', 2); hold on;

  extendPolyLoading(poly_00, color_00, 1.0);
  extendPolyUnloading(poly_09, color_09, 1.0);
  extendPolyUnloading(poly_18, color_18, 1.0);
  extendPolyUnloading(poly_27, color_27, 1.0);
  extendPolyUnloading(poly_36, color_36, 1.0);
  extendPolyUnloading(poly_45, color_45, 1.0);
  

  p = 1.0e8;
  dp = 0.1 * p;
  seg = [[-1 p - dp];[1 p - dp]];
  plot(seg(:,1), seg(:,2), 'r--', 'linewidth', 2);

  [p_00_lo, t_00_lo] = lineIntersectSegLinear(seg, poly_00);
  [p_09_lo, t_09_lo] = lineIntersectSegLinear(seg, poly_09);
  [p_18_lo, t_18_lo] = lineIntersectSegLinear(seg, poly_18);
  [p_27_lo, t_27_lo] = lineIntersectSegLinear(seg, poly_27);
  [p_36_lo, t_36_lo] = lineIntersectSegLinear(seg, poly_36);
  [p_45_lo, t_45_lo] = lineIntersectSegLinear(seg, poly_45);
  plot(p_00_lo(1), p_00_lo(2), 'kx', 'linewidth', 2, 'markersize', 8);
  plot(p_09_lo(1), p_09_lo(2), 'kx', 'linewidth', 2, 'markersize', 8);
  plot(p_18_lo(1), p_18_lo(2), 'kx', 'linewidth', 2, 'markersize', 8);
  plot(p_27_lo(1), p_27_lo(2), 'kx', 'linewidth', 2, 'markersize', 8);
  plot(p_36_lo(1), p_36_lo(2), 'kx', 'linewidth', 2, 'markersize', 8);
  plot(p_45_lo(1), p_45_lo(2), 'kx', 'linewidth', 2, 'markersize', 8);

  seg = [[-1 p + dp];[1 p + dp]];
  [p_00_hi, t_00_hi] = lineIntersectSegLinear(seg, poly_00);
  [p_09_hi, t_09_hi] = lineIntersectSegLinear(seg, poly_09);
  [p_18_hi, t_18_hi] = lineIntersectSegLinear(seg, poly_18);
  [p_27_hi, t_27_hi] = lineIntersectSegLinear(seg, poly_27);
  [p_36_hi, t_36_hi] = lineIntersectSegLinear(seg, poly_36);
  [p_45_hi, t_45_hi] = lineIntersectSegLinear(seg, poly_45);
  plot(seg(:,1), seg(:,2), 'b--', 'linewidth', 2);
  %plot(p_00(1), p_00(2), 'kx', 'linewidth', 2, 'markersize', 8);
  %plot(p_09(1), p_09(2), 'kx', 'linewidth', 2, 'markersize', 8);
  %plot(p_18(1), p_18(2), 'kx', 'linewidth', 2, 'markersize', 8);
  %plot(p_27(1), p_27(2), 'kx', 'linewidth', 2, 'markersize', 8);
  %plot(p_36(1), p_36(2), 'kx', 'linewidth', 2, 'markersize', 8);
  %plot(p_45(1), p_45(2), 'kx', 'linewidth', 2, 'markersize', 8);

  fig2 = figure();
  eps_p = [0, 0.06485046787894522, 0.15676570118384126, 0.2560980039246019, 0.35941360946433987, 0.3856727005415841];
  eps_v_lo = [p_00_lo(1), p_09_lo(1), p_18_lo(1), p_27_lo(1), p_36_lo(1), p_45_lo(1)];
  eps_v_hi = [p_00_hi(1), p_09_hi(1), p_18_hi(1), p_27_hi(1), p_36_hi(1), p_45_hi(1)];
  plot(eps_v_lo, eps_p, '-x', 'color', [1 0 0], 'linewidth', 2); hold on;
  plot(eps_v_hi, eps_p, '-x', 'color', [0 0 1], 'linewidth', 2); hold on;

  seg = [[0.2 -1];[0.2 1]];
  poly_lo = [eps_v_lo' eps_p'];
  poly_hi = [eps_v_hi' eps_p'];
  [eps_lo, t_lo] = lineIntersectSegLinear(seg, poly_lo)
  [eps_hi, t_hi] = lineIntersectSegLinear(seg, poly_hi)
  plot(eps_lo(1), eps_lo(2), 'kx', 'linewidth', 2, 'markersize', 8);
  plot(eps_hi(1), eps_hi(2), 'kx', 'linewidth', 2, 'markersize', 8);

  dp = 2 * dp;
  deps_v = eps_hi(2) - eps_lo(2);
  K = dp/deps_v
  

end

function [poly] = changeUnits(poly)
  poly(:,1) = poly(:,1)*1.0e-2;
  poly(:,2) = poly(:,2)*1.0e6;
end

function extendPolyLoading(poly, color, eps_max)
  [n, m] = size(poly);
  %eps_max = 1.0;
  p_n_1 = poly(n-1, 2);
  p_n = poly(n, 2);
  eps_n_1 = poly(n-1, 1);
  eps_n = poly(n, 1);
  t = (eps_max - eps_n_1)/(eps_n - eps_n_1);
  p_max = (1 - t)*p_n_1 + t*p_n;
  plot([eps_n eps_max],[p_n p_max], '--', 'color', color, 'linewidth', 2); hold on;
end

function extendPolyUnloading(poly, color)
  [n, m] = size(poly);
  eps_max = 1.0;
  p_n_1 = poly(1, 2);
  p_n = poly(2, 2);
  eps_n_1 = poly(1, 1);
  eps_n = poly(2, 1);
  t = (eps_max - eps_n_1)/(eps_n - eps_n_1);
  p_max = (1 - t)*p_n_1 + t*p_n;
  plot([eps_n eps_max],[p_n p_max], '--', 'color', color, 'linewidth', 2); hold on;
end

function nonconvex_poly()
  poly = [0,0,0; 1,2,0; 2,3,0; 3,2,0; 4,3,0; 5,0,0];
  plot(poly(:,1), poly(:,2), 'ro-', 'linewidth', 2); hold on;
  grid on;
  axis equal;

  seg = [-1,1,0; -1,2,0];
  %plotSegAndIntersect(seg, poly);

  seg = [1,1,0; 2,2,0];
  %plotSegAndIntersect(seg, poly);

  seg = [1,1,0; 1,3,0];
  %plotSegAndIntersect(seg, poly);

  seg = [2.5,2,0; 2.5,3,0];
  %plotSegAndIntersect(seg, poly);

  seg = [2,1,0; 0,3,0];
  %plotSegAndIntersect(seg, poly);

  seg = [2,0,0; 1,1,0];
  %plotSegAndIntersect(seg, poly);

  seg = [4.33081187259450e-01, 9.95411511939463e-01, 0;...
         1.37170358414487e+00, 2.44102702831043e+00, 0];
  plotSegAndIntersect(seg, poly);
end

function plotSegAndIntersect(seg, poly)
  plotSeg(seg, [0.1, 1, 0]);
  [p_lo, p_hi, t_lo, t_hi, s1_lo, s1_hi, s2_lo, s2_hi] = lineIntersectSegBinary(seg, poly)
  plot(p_lo(1), p_lo(2), 'kx', 'MarkerSize', 10);
  plot(p_hi(1), p_hi(2), 'mx', 'MarkerSize', 10);
end

function plotSeg(seg, color)
  plot(seg(:,1), seg(:,2), '-', 'color', color, 'linewidth', 2);
end
function plotSegThick(seg, color)
  plot(seg(:,1), seg(:,2), '-x', 'color', color, 'linewidth', 5, 'markersize', 10);
end

function [q, t, s1, s2] = lineIntersectSegAndSide(line, seg)
  p0 = line(1,:)';
  p1 = line(2,:)';
  q0 = seg(1,:)';
  q1 = seg(2,:)';
  p1p0 = p1 - p0;
  q1q0 = q1 - q0;
  J = [p1p0 -q1q0];
  f = (q0 - p0);

  J1 = [(q0 - p0)  (q1 - p0)];
  J2 = [(q0 - p1)  (q1 - p1)];
  s1 = sign(det(J1(1:2,1:2)));
  s2 = sign(det(J2(1:2,1:2)));

  if (det(J(1:2,1:2)) == 0) 
    % Parallel or collinear
    q1q0 = q1 - q0;
    p0q0 = p0 - q0;
    if (p0q0 == 0) 
      p0q0 = p1 - q0;
    endif  
    angle = dot(q1q0/norm(q1q0), p0q0/norm(p0q0))
    if (abs(angle) < 1.0e-6 || abs(abs(angle)-1) < 1.0e-6)
      q0
      q1
      tx0 = (p0(1) - q0(1))/(q1(1) - q0(1));
      ty0 = (p0(2) - q0(2))/(q1(2) - q0(2));
      tx1 = (p1(1) - q0(1))/(q1(1) - q0(1));
      ty1 = (p1(2) - q0(2))/(q1(2) - q0(2));
      txy = [tx0 ty0 tx1 ty1]
      if (tx0 == ty0 && tx1 == ty1)
        if ((tx0 < 0 || tx0 > 1) && (tx1 < 0 || tx1 > 1))
          t = [-1 -1];
        else
          tx0 = max(0, tx0)
          tx1 = max(0, tx1)
          tx0 = min(1, tx0)
          tx1 = min(1, tx1)
          tx = 0.5 * (tx0 + tx1);
          t = [tx tx];
        endif
      else
        t = [-1 -1];
      endif
    else
      t = [-1 -1];
    endif
    q = (1 - t(2)) * q0 + t(2) * q1;
  else
    t =  inv(J(1:2,1:2)) * f(1:2);
    q = (1 - t(2)) * q0 + t(2) * q1;
  endif

end



function [p, t] = lineIntersectSegLinear(seg, poly)

  for i=1:length(poly)-1
    seg_poly = [poly(i,:); poly(i+1,:)];
    %seg;
    %seg_poly;
    [p, t, s1, s2] = lineIntersectSegAndSide(seg, seg_poly);
    if (!(t(2) < 0) && !(t(2) > 1))
      %plot(p(1), p(2), 'kx', 'MarkerSize', 10);
      return;
    else
      if (i == 1 && t(2) < 0)
        plot(p(1), p(2), 'rx', 'MarkerSize', 10);
        return;
      endif
    endif
  end

end

% Assuming convex poly points are arranged clockwise
function [p_lo, p_hi, t_lo, t_hi, s1_lo, s1_hi, s2_lo, s2_hi] = lineIntersectSegBinary(line_seg, poly)

  [n, dim] = size(poly);
  n_mid = round(n/2);

  %final = [n n_mid poly(1,:) poly(n_mid,:) poly(n,:)]
  seg_mid_hi = [poly(n_mid,:); poly(n, :)];
  if (n == 2)
    [p_lo, t_lo, s1_lo, s2_lo] = lineIntersectSegAndSide(line_seg, poly);
    p_hi = p_lo;
    t_hi = t_lo;
    s1_hi = s1_lo;
    s2_hi = s2_lo;
    plotSegThick(poly, [0, 0, 0]);
    return
  endif

  seg_lo_mid = [poly(1,:); poly(n_mid,:)];
  seg_mid_hi = [poly(n_mid,:); poly(n, :)];
  plotSeg(seg_lo_mid, [1, 0, 1]);
  plotSeg(seg_mid_hi, [1, 1, 0]);
  [p_lo, t_lo, s1_lo, s2_lo] = lineIntersectSegAndSide(line_seg, seg_lo_mid);
  [p_hi, t_hi, s1_hi, s2_hi] = lineIntersectSegAndSide(line_seg, seg_mid_hi);
  %plot(p_lo(1), p_lo(2), 'kx', 'MarkerSize', 10);
  %plot(p_hi(1), p_hi(2), 'kx', 'MarkerSize', 10);


  % Recurse on low side
  if ((s1_lo > 0 && s2_lo > 0) || (s1_hi < 0 && s2_hi < 0))
    lo1 = [n n_mid s1_lo s2_lo s1_hi s2_hi]
    [p_lo, p_hi, t_lo, t_hi, s1_lo, s1_hi, s2_lo, s2_hi] = lineIntersectSegBinary(line_seg, poly(1:n_mid,:));
  elseif ((s1_lo < 0 && s2_lo < 0) || (s1_hi > 0 && s2_hi > 0))
    lo2 = [n n_mid s1_lo s2_lo s1_hi s2_hi]
    [p_lo, p_hi, t_lo, t_hi, s1_lo, s1_hi, s2_lo, s2_hi] = lineIntersectSegBinary(line_seg, poly(n_mid:n,:));
  else
    mid = [n n_mid s1_lo s2_lo s1_hi s2_hi]
    if (t_lo(2) < 1)
      [p_lo, p_hi, t_lo, t_hi, s1_lo, s1_hi, s2_lo, s2_hi] = lineIntersectSegBinary(line_seg, poly(1:n_mid,:));
    else
      [p_lo, p_hi, t_lo, t_hi, s1_lo, s1_hi, s2_lo, s2_hi] = lineIntersectSegBinary(line_seg, poly(n_mid:n,:));
    endif
  end


end

