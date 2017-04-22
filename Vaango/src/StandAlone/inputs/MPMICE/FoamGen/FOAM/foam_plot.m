  clear all;
  close all;
  load foam_40mm_v3.1.vf

  x = foam_40mm_v3_1(:,1);
  y = foam_40mm_v3_1(:,2);
  vf = foam_40mm_v3_1(:,4);

  clear foam_40mm_v3_1;

  for i=1:size(x)
    if (vf(i) > 0.2) 
      plot(x(i),y(i),'r.'); hold on;
    end
  end

  load foam_40mm_v3.1.pts
  xpt = foam_40mm_v3_1(:,5);
  ypt = foam_40mm_v3_1(:,6);
  plot(xpt, ypt, 'bx'); hold on;


