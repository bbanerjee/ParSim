function test

  load Force40mm20ms300K_ng.dat
  load Force40mm_ng_200ms.dat

  force20 =  Force40mm20ms300K_ng;
  force200 =  Force40mm_ng_200ms;

  v20 = 20.0;
  v200 = 200.0;
  l = 40.0e-3;
  w = 2.0e-4;

  t20 = force20(:,1);
  fy20 = force20(:,3);
  t200 = force200(:,1);
  fy200 = force200(:,3);

  sig20 = fy20/(l*w);
  sig200 = fy200/(l*w);

  d20 = t20*v20;
  d200 = t200*v200;
  lnew20 = 1.0 - d20/l;
  lnew200 = 1.0 - d200/l;
   
  eps20 = log(lnew20);
  eps200 = log(lnew200);

  figure;
  p1 = plot(-eps20, sig20*1.0e-6, 'b-', 'LineWidth', 2); hold on;
  p2 = plot(-eps200, sig200*1.0e-6, 'r-', 'LineWidth', 2); hold on;
  xlabel('True strain', 'FontName', 'times', 'FontSize', 18);
  ylabel('True Stress (MPa)', 'FontName', 'times', 'FontSize', 18);
  set(gca, 'LineWidth', 3,  'FontName', 'times', 'FontSize', 18);
  legend([p1, p2], 'v = 20 m/s', 'v = 200 m/s');


  load Force40mm20ms300K_ng.dat
  load Force40mm20ms450K_ng.dat

  force300 =  Force40mm20ms300K_ng;
  force450 =  Force40mm20ms450K_ng;

  v = 20.0;
  l = 40.0e-3;
  w = 2.0e-4;

  t300 = force300(:,1);
  fy300 = force300(:,3);
  t450 = force450(:,1);
  fy450 = force450(:,3);

  sig300 = fy300/(l*w);
  sig450 = fy450/(l*w);

  d300 = t300*v;
  d450 = t450*v;
  lnew300 = 1.0 - d300/l;
  lnew450 = 1.0 - d450/l;
   
  eps300 = log(lnew300);
  eps450 = log(lnew450);

  figure;
  p1 = plot(-eps300, sig300*1.0e-6, 'b-', 'LineWidth', 2); hold on;
  p2 = plot(-eps450, sig450*1.0e-6, 'r-', 'LineWidth', 2); hold on;
  xlabel('True strain', 'FontName', 'times', 'FontSize', 18);
  ylabel('True Stress (MPa)', 'FontName', 'times', 'FontSize', 18);
  set(gca, 'LineWidth', 3,  'FontName', 'times', 'FontSize', 18);
  title('Impact Velocity = 20 m/s (~ 500/s)', ...
         'FontName', 'times', 'FontSize', 18);
  legend([p1, p2], 'T = 300 K', 'T = 450 K');
