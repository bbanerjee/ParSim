  clear all;
  close all;
  load foam_v4.2.vf

  x = foam_v4(:,1);
  y = foam_v4(:,2);
  vf = foam_v4(:,4);

  size(x)
  for i=1:size(x)
    if (vf(i) > 1.0e-5) 
      plot(x(i),y(i),'r.'); hold on;
    end
  end

