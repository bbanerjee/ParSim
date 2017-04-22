function part_plot
  load part_20mm.dat
  load part_40mm.dat
  load part_60mm.dat

  plotAll(part_20mm, 20);
  plotAll(part_40mm, 40);
  plotAll(part_60mm, 60);

function plotAll(part, size)

  figure;

  r = part(:,2);
  xc = part(:,5);
  yc = part(:,6);
  for i=1:length(r)
    plotCircle(r(i),xc(i),yc(i));
  end
  set(gca, 'XLim', [0 size], 'YLim', [0 size]);
  axis square;
  set(gca, 'LineWidth', 3, 'FontName', 'times', 'FontSize', 18);
  xlabel('mm',  'FontName', 'times', 'FontSize', 18);
  ylabel('mm',  'FontName', 'times', 'FontSize', 18);
  

function plotCircle(r,xc,yc)

  thetamin = 0;
  thetamax = 360.0;
  ntheta = 64;
  dtheta = (thetamax-thetamin)/ntheta;
  for i=1:ntheta+1
    theta = (i-1)*dtheta;
    x(i) = xc + r*cos(theta);
    y(i) = yc + r*sin(theta);
  end
  plot(x*1000, y*1000, 'b-', 'LineWidth', 2); hold on;

