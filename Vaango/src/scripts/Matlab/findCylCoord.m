function findCylCoord

  %pos = [0 0 2]
  %bot = [1 1 0]
  %top = [1.5 1.5 1]
  %findCoord(pos, top, bot)

  %pos = [2 1 0.5]
  %bot = [1 1 0]
  %top = [1 1 1]
  %findCoord(pos, top, bot)

  %r = 1;
  %theta = 30;
  %c = r*cos(theta*pi/180);
  %s = r*sin(theta*pi/180);
  %pos = [1+c 1+s 0.5]
  %bot = [1 1 0]
  %top = [1 1 1]
  %findCoord(pos, top, bot)
  
  %r = 2;
  %theta = 60;
  %c = r*cos(theta*pi/180);
  %s = r*sin(theta*pi/180);
  %pos = [1+c 1+s 0.5]
  %bot = [1 1 1]
  %top = [1 1 2]
  %findCoord(pos, top, bot)

  %r = 2;
  %theta = 0;
  %c = r*cos(theta*pi/180);
  %s = r*sin(theta*pi/180);
  %pos = [1+c 1+s 0.5]
  %bot = [1 1 1]
  %top = [1 1 2]
  %tilt_angle = -30*pi/180
  %rot_axis = [0 1 0]
  %[tilted_axis_top] = findTiltedAxis(bot, top, tilt_angle, rot_axis)
  %[rot_point] = rotateLocal(pos - bot, tilted_axis_top - bot)
  %rot_point = rot_point + bot
  %findCoord(rot_point, tilted_axis_top, bot)
  
  %right = [2 1 1];
  %[rot_right] = rotateLocal(right - bot, tilted_axis_top - bot)
  %rot_right = rot_right + bot

  %figure
  %plot3([bot(1) top(1)], [bot(2) top(2)], [bot(3) top(3)], 'r-', 'LineWidth', 2); hold on
  %grid on;
  %plot3(0, 0, 0, 'bo', 'MarkerSize', 10);
  %plot3(pos(1), pos(2), pos(3), 'rx', 'MarkerSize', 10);
  %plot3(right(1), right(2), right(3), 'ro', 'MarkerSize', 10);
  %plot3([bot(1) tilted_axis_top(1)], [bot(2) tilted_axis_top(2)], [bot(3) tilted_axis_top(3)], 'm-', 'LineWidth', 2); hold on
  %plot3(rot_point(1), rot_point(2), rot_point(3), 'mx', 'MarkerSize', 10);
  %plot3(rot_right(1), rot_right(2), rot_right(3), 'mo', 'MarkerSize', 10);
  %xlabel('1');
  %ylabel('2');
  %zlabel('3');

  %r = 2;
  %theta = 125;
  %c = r*cos(theta*pi/180);
  %s = r*sin(theta*pi/180);
  %pos = [1.5+c 0.5+s 0.5]
  %bot = [1.5 0.5 1]
  %top = [1.5 0.5 2]
  %tilt_angle = -70*pi/180
  %rot_axis = [0 1 0]
  %[tilted_axis_top] = findTiltedAxis(bot, top, tilt_angle, rot_axis)
  %[rot_point] = rotateLocal(pos - bot, tilted_axis_top - bot)
  %rot_point = rot_point + bot
  %findCoord(rot_point, tilted_axis_top, bot)
  
  %right = [2 1 1];
  %[rot_right] = rotateLocal(right - bot, tilted_axis_top - bot)
  %rot_right = rot_right + bot

  %figure
  %plot3([bot(1) top(1)], [bot(2) top(2)], [bot(3) top(3)], 'r-', 'LineWidth', 2); hold on
  %grid on;
  %plot3(0, 0, 0, 'bo', 'MarkerSize', 10);
  %plot3(pos(1), pos(2), pos(3), 'rx', 'MarkerSize', 10);
  %plot3(right(1), right(2), right(3), 'ro', 'MarkerSize', 10);
  %plot3([bot(1) tilted_axis_top(1)], [bot(2) tilted_axis_top(2)], [bot(3) tilted_axis_top(3)], 'm-', 'LineWidth', 2); hold on
  %plot3(rot_point(1), rot_point(2), rot_point(3), 'mx', 'MarkerSize', 10);
  %plot3(rot_right(1), rot_right(2), rot_right(3), 'mo', 'MarkerSize', 10);
  %xlabel('1');
  %ylabel('2');
  %zlabel('3');

  r = 2;
  theta = 30;
  c = r*cos(theta*pi/180);
  s = r*sin(theta*pi/180);
  pos = [1+c 1+s 1.5]
  bot = [1 1 1]
  top = [1 1 2]
  right = [4 1 1];
  
  
  figure
  %[xx yy zz] = cylinder;
  %hh = surf(bot(1)+xx*r, bot(2)+yy*r, bot(3)+zz*(top(3)-bot(3))); hold on;
  plot3([bot(1) top(1)], [bot(2) top(2)], [bot(3) top(3)], 'r-', 'LineWidth', 2); hold on
  plot3(0, 0, 0, 'bo', 'MarkerSize', 10);
  plot3(pos(1), pos(2), pos(3), 'rx', 'MarkerSize', 10);
  plot3([bot(1) right(1)], [bot(3) right(2)], [bot(3) right(3)], 'r--', 'MarkerSize', 10);
  xlabel('1');
  ylabel('2');
  zlabel('3');
  grid on;

  tilt_angle = -60*pi/180
  rot_axis = [0 1 0]
  [tilted_axis_top] = findTiltedAxis(bot, top, tilt_angle, rot_axis)
  [rot_point] = rotateLocal(pos - bot, tilted_axis_top - bot)
  rot_point = rot_point + bot
  [rot_right] = rotateLocal(right - bot, tilted_axis_top - bot)
  rot_right = rot_right + bot

  plot3([bot(1) tilted_axis_top(1)], [bot(2) tilted_axis_top(2)], [bot(3) tilted_axis_top(3)], 'm-', 'LineWidth', 2); hold on
  plot3(rot_point(1), rot_point(2), rot_point(3), 'mx', 'MarkerSize', 10);
  plot3([bot(1) rot_right(1)], [bot(1) rot_right(2)], [bot(1) rot_right(3)], 'm--', 'MarkerSize', 10);
  %[xx yy zz] = cylinder;
  %hh1 = mesh(bot(1)+xx*r, bot(2)+yy*r, bot(3)+zz*(top(3)-bot(3))); hold on;
  %rotate(hh1, [0 1 0], 60, bot);
  axis equal

  tilt_angle = 45*pi/180
  rot_axis = [0 0 1]
  [tilted_axis_top] = findTiltedAxis(bot, tilted_axis_top, tilt_angle, rot_axis);
  [rot_point] = rotateLocal(pos - bot, tilted_axis_top - bot)
  rot_point = rot_point + bot
  [rot_right] = rotateLocal(right - bot, tilted_axis_top - bot)
  rot_right = rot_right + bot

  plot3([bot(1) tilted_axis_top(1)], [bot(2) tilted_axis_top(2)], [bot(3) tilted_axis_top(3)], 'g-', 'LineWidth', 2); hold on
  plot3(rot_point(1), rot_point(2), rot_point(3), 'gx', 'MarkerSize', 10);
  plot3([bot(1) rot_right(1)], [bot(2) rot_right(2)], [bot(3) rot_right(3)], 'g--', 'MarkerSize', 10);

  zaxis = [0 0 1]
  rangle = acos(dot(tilted_axis_top - bot, zaxis)/norm(tilted_axis_top - bot))
  raxis = cross(tilted_axis_top - bot, zaxis)
  
  [xx yy zz] = cylinder;
  hh2 = surf(bot(1)+xx*r, bot(2)+yy*r, bot(3)+zz*(top(3)-bot(3))); hold on;
  rotate(hh2, raxis, -rangle*180/pi, bot);

  findCoord(rot_point, tilted_axis_top, bot)
  


function findCoord(pos, top, bot)

  axis_e3 = [0 0 1]
  axis_ez = top - bot;
  axis_ez = axis_ez/norm(axis_ez)

  nn = axis_ez'*axis_ez
  One = [[1 0 0];[0 1 0];[0 0 1]]

  particleLoc = pos - bot
  particleProj = nn*particleLoc'
  %zz = norm(particleProj)
  
  bot_shift = bot - bot;
  top_shift = top - bot;
  top_shift - bot_shift
  [max_val, max_index] = max(top_shift - bot_shift)
  t = (particleProj(max_index) - bot_shift(max_index))/(top_shift(max_index) - bot_shift(max_index))
  zmax = norm(top_shift)
  zz = t*zmax

  projMatrix = One - nn
  particleProj = projMatrix*particleLoc'
  r = norm(particleProj)

  axis_e1 = [1  0  0];
  axisLoc = axis_e1
  axisProj = projMatrix*axisLoc'

  rvec = particleProj/r
  thetavec = axisProj/norm(axisProj)

  theta = acos(dot(rvec, thetavec))
  theta = theta*180/pi

  rvec = 4*rvec + bot';
  thetavec = 4*thetavec + bot';
  plot3([bot(1) rvec(1)],[bot(2) rvec(2)],[bot(3) rvec(3)], 'LineWidth', 3, 'Color', [255/255 153/255 51/255]);
  plot3([bot(1) thetavec(1)],[bot(2) thetavec(2)],[bot(3) thetavec(3)], 'LineWidth', 3, 'Color', [51/255 153/255 255/255]);

function [tilted_axis_point] = findTiltedAxis(bot, top, angle, rot_axis)

  point = top - bot
  rot_angle = angle
  %rot_axis = [0 1 0]
  R = rotation_matrix(rot_angle, rot_axis)
  tilted_axis_point = (R*point')'
  tilted_axis_point = tilted_axis_point + bot 
  
function [rot_point] = rotateLocal(point, tilted_axis)

  zaxis = [0 0 1]
  rot_angle = acos(dot(tilted_axis, zaxis)/norm(tilted_axis))
  rot_axis = cross(tilted_axis, zaxis)

  R = rotation_matrix(rot_angle, rot_axis)
  rot_point = (R*point')'

function [R] = rotation_matrix(angle, axis)

   axis = axis/norm(axis);
   ca = cos(angle);
   sa = sin(angle);
   I = [[1 0 0];[0 1 0];[0 0 1]];
   aa = axis'*axis;
   A = [[0 axis(3) -axis(2)];[-axis(3) 0 axis(1)];[axis(2) -axis(1) 0]];
   R = (I - aa)*ca + aa + A*sa;
