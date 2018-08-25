function table_intersect_spline()

  format long e

  z = [0, 1, 2, 3, 4, 5, 6, 7];
  rprime = [0, 2, 3, 2, 3, 0, 0, 0];

  plot(z, rprime, 'r-', 'LineWidth', 2); hold on;
  axis equal
  grid on
  %axis square

  %test1(z, rprime);
  %test2(z, rprime);
  %test3(z, rprime);
  %test4(z, rprime);
  test5(z, rprime);
end

function test1(z, rprime)

  seg_z = [1 2];
  seg_rprime = [1 2];
  int_line = [0.0857864 0.0857864];
  int_spline = [0 0];

  plot(seg_z, seg_rprime, 'bx-', 'LineWidth', 3, 'Markersize', 7);
  plot(int_line(1), int_line(2), 'gs', 'LineWidth', 1, 'Markersize', 10);
  plot(int_spline(1), int_spline(2), 'ms', 'LineWidth', 1, 'Markersize', 10);
  plotBSpline(z, rprime);

end

function test2(z, rprime)

  seg_z = [2 0];
  seg_rprime = [1 3];
  int_line = [1 2];
  int_spline = [1.05051 1.94949];

  plot(seg_z, seg_rprime, 'bx-', 'LineWidth', 3, 'Markersize', 7);
  plot(int_line(1), int_line(2), 'gs', 'LineWidth', 1, 'Markersize', 10);
  plot(int_spline(1), int_spline(2), 'ms', 'LineWidth', 1, 'Markersize', 10);
  plotBSpline(z, rprime);

end

function test3(z, rprime)

  seg_z = [2 1];
  seg_rprime = [0 1];
  int_line = [6.5 -4.5];
  int_spline = [2 0];

  plot(seg_z, seg_rprime, 'bx-', 'LineWidth', 3, 'Markersize', 7);
  plot(int_line(1), int_line(2), 'gs', 'LineWidth', 1, 'Markersize', 10);
  plot(int_spline(1), int_spline(2), 'ms', 'LineWidth', 1, 'Markersize', 10);
  plotBSpline(z, rprime);

end

function test4(z, rprime)
  seg_z = [4.33081187259450e-01 1.37170358414487e+00];
  seg_rprime = [9.95411511939463e-01 2.44102702831043e+00];
  int_line = [0.714147 1.42829];
  int_spline = [1.08032 1.99226];

  plot(seg_z, seg_rprime, 'bx-', 'LineWidth', 3, 'Markersize', 7);
  plot(int_line(1), int_line(2), 'gs', 'LineWidth', 1, 'Markersize', 10);
  plot(int_spline(1), int_spline(2), 'ms', 'LineWidth', 1, 'Markersize', 10);
  plotBSpline(z, rprime);

end

function test5(z, rprime)
  seg_z = [1.5 2.5];
  seg_rprime = [2.5 2.5];
  int_line = [1.5 2.5];
  int_spline = [1.5 2.5];

  plot(seg_z, seg_rprime, 'bx-', 'LineWidth', 3, 'Markersize', 7);
  plot(int_line(1), int_line(2), 'gs', 'LineWidth', 1, 'Markersize', 10);
  plot(int_spline(1), int_spline(2), 'ms', 'LineWidth', 1, 'Markersize', 10);
  plotBSpline(z, rprime);

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

function plotBSpline(z, rprime)

  for index = 2:length(z)-1
    loc = 'mid';
    p1 = [z(index-1) rprime(index-1)];
    p2 = [z(index) rprime(index)];
    p3 = [z(index+1) rprime(index+1)];
    t = linspace(0, 1, 20);
    for i=1:length(t)
      [B, T, N] = computeBSpline(loc, p1, p2, p3, t(i));
      plot(B(1), B(2), 'k.');
    end
  end
end

