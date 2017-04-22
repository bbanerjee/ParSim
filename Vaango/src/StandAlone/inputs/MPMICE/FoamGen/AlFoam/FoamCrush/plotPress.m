function plotPress

  load pressCC_ng.dat;

  l = 0.04;
  w = 2.0e-4;
  nc = 201;
  dl = l/nc;
  ac = dl*w;
  press = pressCC_ng(:,4);
  time = pressCC_ng(:,1);

  press = press - 101325.0;

  fy = press*ac;

  count = 1;
  fysum(count) = 0;
  for i=1:length(fy)
    fysum(count) = fysum(count) + fy(i);
    if (time(i) == 200) 
      count = count+1;
      fysum(count) = 0;
    end
  end
  length(fysum)

  for i=1:length(fysum)
    sig(i) = fysum(i)/(l*w);
  end
  length(sig)

  load timev1.dat
  t = timev1(:,2);

  v = 200.0;

  d = t*v;
  lnew = l - d;
   
  for i=1:length(lnew)
    eps(i) = log(lnew(i)/l);
  end
  length(eps)

  plot(-eps(1:length(eps)-2), sig(1:length(eps)-2)*1.0e-5, 'b-', 'LineWidth', 2);
  xlabel('True strain', 'FontName', 'times', 'FontSize', 18);
  ylabel('True Stress (MPa)', 'FontName', 'times', 'FontSize', 18);
  set(gca, 'LineWidth', 3,  'FontName', 'times', 'FontSize', 18);
