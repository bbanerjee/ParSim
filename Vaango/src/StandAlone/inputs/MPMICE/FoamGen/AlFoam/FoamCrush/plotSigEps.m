function plotSigEps

  load Force40mm_g_200ms.dat
  load Force40mm_ng_200ms.dat

  force_g =  Force40mm_g_200ms;
  force_ng =  Force40mm_ng_200ms;

  v = 200.0;
  l = 40.0e-3;
  w = 2.0e-4;

  t_g = force_g(:,1);
  fy_g = force_g(:,3);
  t_ng = force_ng(:,1);
  fy_ng = force_ng(:,3);

  sig_g = fy_g/(l*w);
  sig_ng = fy_ng/(l*w);

  d_g = t_g*v;
  d_ng = t_ng*v;
  lnew_g = l - d_g;
  lnew_ng = l - d_ng;
   
  for i=1:length(lnew_g)
     eps_g(i) = log(lnew_g(i)/l);
  end
  for i=1:length(lnew_ng)
     eps_ng(i) = log(lnew_ng(i)/l);
  end

  figure;
  p1 = plot(-eps_ng, sig_ng*1.0e-6, 'r-', 'LineWidth', 2); hold on;
  p2 = plot(-eps_g, sig_g*1.0e-6, 'b-', 'LineWidth', 2);
  xlabel('True strain', 'FontName', 'times', 'FontSize', 18);
  ylabel('True Stress (MPa)', 'FontName', 'times', 'FontSize', 18);
  set(gca, 'LineWidth', 3,  'FontName', 'times', 'FontSize', 18);
  title('Impact Velocity = 200 m/s (~ 5000/s)', ...
         'FontName', 'times', 'FontSize', 18);
  legend([p1, p2], 'No Gas', 'With Gas');
