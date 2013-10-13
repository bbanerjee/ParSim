function test_ordering

  close all;
  clear all;

  load bunny1.txt
  x = bunny1(:,1);
  y = bunny1(:,2);
  plot(x, y, 'bx'); hold on;
  plot(x, y, 'r-'); hold on;

  angleApproach(x,y)

function angleApproach(x,y)

  % Find mean x and y
  x_mean = mean(x);
  y_mean = mean(y);

  % Draw lines from mean to all points
  % for ii=1:length(x)
  %   plot([x_mean x(ii)],[y_mean y(ii)],'b-'); 
  % end

  % Take an input point
  [x_inp, y_inp] = ginput(1);

  % Draw lines from input pt to all points
  for ii=1:length(x)
    plot([x_inp x(ii)],[y_inp y(ii)],'b-'); 
    ang(ii,1) = atan2(y(ii)-y_inp, x(ii)-x_inp);
  end

  % place in matrix and sort
  mat = [x y ang]
  mat = sortrows(mat, 3);
  x_out = mat(:,1);
  y_out = mat(:,2);
  hold off; 
  plot(x, y, 'bx'); hold on;
  plot(x_out, y_out, 'm-'); hold on;

  % Copy back
  x = x_out;
  y = y_out;

  x(length(x)+1) = x(1);
  y(length(y)+1) = y(1);
  small = true;
  while (small)
    forward = true;
    [xnew, ynew, small] = removeSmallAngles(x,y, forward);
    x = xnew;
    y = ynew;
    plot(x, y, 'b-'); hold on;
    plot(x, y, 'ro'); hold on;
    pause;
    reverse = false;
    [xnew, ynew, small] = removeSmallAngles(x,y, reverse);
    x = xnew;
    y = ynew;
    hold off; plot(x, y, 'b-'); hold on;
    plot(x, y, 'ro'); hold on;
    pause;
  end

  % Save data
  data = [x; y];
  fid = fopen('bunny_sorted.dat','w');
  for ii=1:length(x)
    fprintf(fid,'%12.8f %12.8f\n', x(ii), y(ii));
  end
  fclose(fid);

function [xnew, ynew, small] = removeSmallAngles(x,y, forward)

  small = false;
  xnew = x;
  ynew = y;
  [length(x)]

  % create line segments
  if (forward)
    for ii=1:length(x)-1
      line(ii).x1 = x(ii);
      line(ii).x2 = x(ii+1);
      line(ii).y1 = y(ii);
      line(ii).y2 = y(ii+1);
      line(ii).a_index = ii;
      line(ii).b_index = ii+1;
    end
  else
    count = 0;
    for ii=length(x):-1:2
      count = count+1;
      line(count).x1 = x(ii);
      line(count).x2 = x(ii-1);
      line(count).y1 = y(ii);
      line(count).y2 = y(ii-1);
      line(count).a_index = ii;
      line(count).b_index = ii-1;
    end
  end

  % find angles between adjacent line segments
  for jj=1:length(line)-1
    line1 = line(jj);
    line2 = line(jj+1);
    a_index = line2.a_index;
    b_index = line2.b_index;
    angle = findAngle(line1, line2);
    str = sprintf('%5.1f', angle*180/pi);
    text(x(a_index),y(a_index),str);
    if (angle > 160.0*pi/180.0)
      small = true;
      [angle a_index b_index]
      tempx = xnew(a_index);
      tempy = ynew(a_index);
      xnew(a_index) = xnew(b_index);
      ynew(a_index) = ynew(b_index);
      xnew(b_index) = tempx;
      ynew(b_index) = tempy;
      return;
    end
  end
  

function angle = findAngle(line1, line2)

   a = [line1.x2-line1.x1 line1.y2-line1.y1];
   b = [line2.x2-line2.x1 line2.y2-line2.y1];
   adotb = a(1)*b(1) + a(2)*b(2);
   alen = sqrt(a(1)*a(1)+a(2)*a(2));
   blen = sqrt(b(1)*b(1)+b(2)*b(2));
   angle = acos(adotb/(alen*blen));

function intersectionApproach(x,y)

  x(length(x)+1) = x(1);
  y(length(y)+1) = y(1);

  intersects = true;
  while (intersects)
    [intersects, xnew, ynew] = checkIntersections(x,y);
    plot(xnew, ynew, 'kx'); hold on;
    plot(xnew, ynew, 'k-'); hold on;
    x = xnew;
    y = ynew;
  end

function [inter, xnew, ynew] = checkIntersections(x,y)

  % copy points
  for ii=1:length(x)
    xnew(ii,1) = x(ii);
    ynew(ii,1) = y(ii);
  end

  % create line segments
  for ii=1:length(x)-1
    line(ii).pt(1,:) = [x(ii) y(ii)]; 
    line(ii).pt(2,:) = [x(ii+1) y(ii+1)]; 
    line(ii).a_index = ii;
    line(ii).b_index = ii+1;
  end

  % check intersections
  inter = false;
  for jj=1:length(line)-1
    line1 = line(jj);
    b_index = line1.b_index;
    for kk=jj+1:length(line)
      line2 = line(kk);
      c_index = line2.a_index;
      intersects = intersect(line1, line2);
      if (intersects)
        inter = true;
        [b_index c_index]
        % switch points
        xnew(b_index,1) = x(c_index);
        ynew(b_index,1) = y(c_index);
        xnew(c_index,1) = x(b_index);
        ynew(c_index,1) = y(b_index);
        [x y xnew ynew]
        return;
      end
    end
  end
  

function [intersects] = intersect(line1, line2)

  intersects = false;

  A = line1.pt(1,:);
  B = line1.pt(2,:);
  C = line2.pt(1,:);
  D = line2.pt(2,:);

  CmP = [C(1)-A(1)  C(2)-A(2)];
  r = [B(1)-A(1)  B(2)-A(2)];
  s = [D(1)-C(1)  D(2)-C(2)];

  CmPxr = CmP(1)*r(2) - CmP(2)*r(1);
  CmPxs = CmP(1)*s(2) - CmP(2)*s(1);
  rxs = r(1)*s(2) - r(2)*s(1);

  if (rxs == 0) 
    return
  end

  rxsr = 1.0/rxs;
  t = CmPxs*rxsr;
  u = CmPxr*rxsr;
 
  if (t > 0.0) && (t < 1.0) && (u > 0.0) && (u < 1.0)
    intersects = true;
  end  

