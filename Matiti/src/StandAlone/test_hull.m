function test_hull

  close all;
  clear all;

  load bunny1.txt
  x = bunny1(:,1);
  y = bunny1(:,2);
  plot(x, y, 'bx'); hold on;

  [xmin, min_index] = min(x);
  [xmax, max_index] = max(x);

  hull_index = convhull(x,y);
  plot(x(hull_index), y(hull_index)); hold on;
  [lower_x, lower_y] = findLowerHull(x, y, min_index, max_index, hull_index);


function [lower_x, lower_y] = findLowerHull(x, y, min_index, max_index, hull_index)

  lower_index = zeros(length(x),1);
  count = 0;
  for ii=1:length(hull_index)
    if (hull_index(ii) >= min_index)
      if (hull_index(ii) <= max_index)
        count = count + 1;
        lower_index(count,1) = hull_index(ii); 
        %[hull_index(ii) min_index max_index]
      end
    end
  end
  count = count+1;
  lower_index(count,1) = lower_index(1,1);
  plot(x(lower_index(1:count)), y(lower_index(1:count)), 'r-'); hold on;
  %lower_index
  lower = zeros(count,2);
  for ii=1:count
    lower(ii,1) = x(lower_index(ii));
    lower(ii,2) = y(lower_index(ii));
  end
  lower
  lower = sort(lower, 2);
  lower
  lower_x = lower(:,1)';
  lower_y = lower(:,2)';
  plot(x(lower_index(1:count)), y(lower_index(1:count)), 'ro'); hold on;
  plot(lower_x, lower_y, 'r-');
  
  
