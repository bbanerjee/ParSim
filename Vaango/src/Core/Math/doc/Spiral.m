function Spiral

  clear all
  close all

  test2DSpiral(1.0, 5);
  test3DSpiral(1.0, 4);
  %test3DSpacedSpiral(1.0, 4);
  test3DSpiralPoints(1.0, 4);
  test3DEqualPartition(1.0, 4);
  test3DEqualPartition(2.0, 8);
  test3DEqualPartition(0.5, 2);

function test2DSpiral(radius, num_radial_layers)

  R = radius/num_radial_layers;

  phi_inc = 2*pi/30;
  phi_max = 3*pi^2*radius/(2.0*R);

  ii_max = phi_max/phi_inc;
  ii = 0:ii_max;

  phi = ii*phi_inc;

  x = radius*phi.*cos(phi);
  y = radius*phi.*sin(phi);

  figure;
  plot(x,y); hold on;
  plot(x,y,'r.'); hold on;
  axis equal;
  grid on;

function test3DSpiral(radius, num_radial_layers)

  R = radius/num_radial_layers;

  phi_inc = 2*pi/30;
  phi_max = 3*pi^2*radius/(2.0*R);

  ii_max = phi_max/phi_inc;
  ii = 0:ii_max;

  phi = ii*phi_inc;

  x = radius*cos(phi).*cos(-pi/2 + phi*(pi/phi_max));
  y = radius*sin(phi).*cos(-pi/2 + phi*(pi/phi_max));
  z = -radius*sin(-pi/2 + phi*(pi/phi_max));

  figure;
  plot3(x,y,z); hold on;
  plot3(x,y,z, 'r.'); hold on;
  axis equal;
  grid on;

  %[xx, yy, zz] = sphere;
  %surf(xx, yy, zz);

function test3DSpacedSpiral(radius, num_radial_layers)

  % Compute the resolution
  R = radius/num_radial_layers;

  % Compute max phi and set up array of indexes
  phi_inc = 2*pi/30;
  phi_max = 3*pi^2*radius/(2.0*R);
  ii_max = ceil(phi_max/phi_inc);
  phi_inc = phi_max/ii_max;
  ii = 0:ii_max;
  phi_0 = ii*phi_inc;

  % Compute spiral length
  mm = -phi_max^2/pi^2;
  [~, s_max] = elliptic123(pi, mm);
  [~, ss_0] = elliptic123(phi_0*pi/phi_max, mm);
  s_max = s_max*radius;
  ss_0 = ss_0*radius;

  % Compute m
  %% mm = phi_max^2/(phi_max^2+pi^2);
  [phi_max max(phi_0) mm]

  % Compute s_tilde array
  ss = s_max*(ii/ii_max)/radius;
  %% s_bar = ss/s_max;
  %% s_tilde = pi^2/sqrt(phi_max^2+pi^2)*s_bar;

  % Invert the elliptic integral to get beta
  format short e
  %beta = inverselliptic2(s_tilde, mm, eps);
  %% [beta_bisect, beta_newton] = ellipticEinv(s_tilde, mm, 0.0, pi/2.0, eps);
  [beta_bisect, beta_newton] = ellipticEinv(phi_max, ss, mm, 0.0, phi_max, eps);
  beta = beta_bisect;
  %% alpha = 2.0*beta;
  %% phi = (alpha - pi/2)*(phi_max/pi);
  phi = beta*phi_max/pi;
  %% [s_tilde' beta' beta_newton' phi' phi_0' ss' ss_0']
  [~, ss_beta] = elliptic123(beta, mm);
  ss_beta = radius*ss_beta;
  [beta' phi' phi_0' ss_beta' ss' ss_0']
  
  x = radius*cos(phi).*cos(-pi/2 + phi*(pi/phi_max));
  y = radius*sin(phi).*cos(-pi/2 + phi*(pi/phi_max));
  z = -radius*sin(-pi/2 + phi*(pi/phi_max));
  x_0 = radius*cos(phi_0).*cos(-pi/2 + phi_0*(pi/phi_max));
  y_0 = radius*sin(phi_0).*cos(-pi/2 + phi_0*(pi/phi_max));
  z_0 = -radius*sin(-pi/2 + phi_0*(pi/phi_max));

  figure;
  plot3(x,y,z, 'r.'); hold on;
  plot3(x,y,z, 'b-'); hold on;
  axis equal;
  grid on;

function test3DSpiralPoints(radius, num_radial_layers)

  R = radius/num_radial_layers;

  phi_inc = 2*pi/30;
  phi_max = 3*pi^2*radius/(2.0*R);

  ii_max = ceil(phi_max/phi_inc);

  center = [0.0;0.0;0.0];

  p = zeros(3, ii_max );

  for ii = 1 : ii_max

    cosphi = ( ( ii_max - ii     ) * ( -1.0 )   ...
             + (     ii - 1 ) * ( +1.0 ) ) ...
             / ( ii_max     - 1 );

    sinphi = sqrt ( 1.0 - cosphi^2 );

    if ( ii == 1 | ii == ii_max )
      theta = 0.0;
    else
      theta = theta + 3.6 / ( sinphi * sqrt ( ii_max ) );
      theta = mod ( theta, 2.0 * pi );
    end

    p(1,ii) = center(1,1) + radius * sinphi * cos ( theta );
    p(2,ii) = center(2,1) + radius * sinphi * sin ( theta );
    p(3,ii) = center(3,1) + radius * cosphi;

  end

  figure;
  x = p(1,:);
  y = p(2,:);
  z = p(3,:);
  plot3(x,y,z, 'r.'); hold on;
  plot3(x,y,z, 'b-'); hold on;
  axis equal;
  grid on;

function test3DEqualPartition(radius, num_radial_layers)

  R = radius/num_radial_layers;

  phi_inc = 2*pi/30;

  % Compute spiral length
  phi_max = 3*pi^2*radius/(2.0*R);
  mm = -phi_max^2/pi^2;
  [~, s_max] = elliptic123(pi, mm);
  s_max = s_max*radius;

  ii_max = ceil(s_max/R);

  points = eq_point_set(2, ii_max);

  figure;
  x = radius*points(1,:);
  y = radius*points(2,:);
  z = radius*points(3,:);
  plot3(x,y,z, 'r.'); hold on;
  plot3(x,y,z, 'b-'); hold on;
  axis equal;
  grid on;
