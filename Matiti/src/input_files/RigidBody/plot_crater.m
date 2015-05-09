function plot_crater

  % Load the data
  dat = load("final_positions");

  % Get the x,y positions
  x = dat(:,1);
  y = dat(:,2);

  % Get the mass and volume
  mass = dat(:,3);
  volume = dat(:,4);

  % Get the limits
  xmin = -0.61;
  xmax = 0.61;
  ymin = -0.46;
  ymax = 0.46;

  % Set grid spacing
  dx = (xmax - xmin)/50;
  nx = floor((xmax - xmin)/dx)+1;
  ny = floor((ymax - ymin)/dx)+1;
  nx
  ny

  % Create a grid 
  xgrid = linspace(xmin, xmax, nx);
  ygrid = linspace(ymin, ymax, ny);

  % Create a Delaunay triangulation
  %tri = delaunay(x, y);

  % Plot the triangulation
  %figure
  %[nn, mm ] = size(tri)
  %for i=1:nn
  %  vx = [x(tri(i,1)) x(tri(i,2)) x(tri(i,3)) x(tri(i,1))];
  %  vy = [y(tri(i,1)) y(tri(i,2)) y(tri(i,3)) y(tri(i,1))];
    %i  
    %vx
  %  plot(vx, vy, "b"); hold on;
  %end
  %plot(x, y, 'ro');
  %axis equal

  %figure
  %triplot(tri, x, y)
  %axis equal

  gmass = zeros(nx,ny);
  gvolume = zeros(nx,ny);
  % Interpolate masses to grid
  for i=1:length(x)
    if (y(i) > ymin) 
      pt = [x(i) y(i)];
      [nodes, weights] = findNodesAndWeights(pt, xmin, ymin, dx);
      for k=1:4
        node = nodes(k,:);
        weight = weights(k);
        gmass(node(1),node(2)) = gmass(node(1),node(2)) + mass(i)*weight;
        gvolume(node(1),node(2)) = gvolume(node(1),node(2)) + volume(i)*weight;
      end
    end
  end

  gdensity = gmass ./ gvolume;
  nan_indices = isnan(gdensity);
  mean_density = mean(mean(gdensity(~nan_indices)));
  gdensity(nan_indices) = mean_density;

  [gx, gy] = meshgrid(xgrid, ygrid);
  %size(gx)
  %size(gy)
  %size(gdensity)
  %max(gdensity)
  %mean(gdensity)
  %min(gdensity)
  figure
  %surf(gx, gy, gdensity');
  surf(gx, gy, gmass');
  shading("interp");

function [nodes, weights] = findNodesAndWeights(pt, xmin, ymin, dx)

  pt_pos = [pt(1) - xmin pt(2) - ymin];
  cell_pos = pt_pos/dx;
  idx = floor(cell_pos);
  nodes = [[idx(1) idx(2)];
           [idx(1)+1 idx(2)];
           [idx(1)+1 idx(2)+1];
           [idx(1) idx(2)+1]];
  NN = cell_pos - idx;
  NN1 = 1.0 - NN;
  weights = [NN1(1)*NN1(2) NN(1)*NN1(2) NN(1)*NN(2) NN1(1)*NN(2)];
  
function [inside] = insideTriangle(pt, v1, v2, v3)

  vec12 = v1 - v2;
  vec32 = v3 - v2;
  vecp2 = pt - v2;
  vec31 = v3 - v1;
  vecp1 = pt - v1;
  vec23 = -vec32;
  vecp3 = pt - v3;

  n = cross(vec12, vec32);
  n2 = cross(vec12, vecp2);
  n1 = cross(vec31, vecp1);
  n3 = cross(vec23, vecp3);

  nn2 = dot(n, n2);
  nn1 = dot(n, n1);
  nn3 = dot(n, n3);

  if (nn1 > 0.0 && nn2 > 0.0 && nn3 > 0.0) 
    inside = true;
  else
    inside = false;
  end

function [value] = interpolateTriangle(pt, v1, v2, v3, q1, q2, q3)

  A = [[v1(1) v1(2) 1];[v2(1) v2(2) 1];[v3(1) v3(2) 1]];
  b = [q1;q2;q3];
  x = inverse(A)*b;

  value = x(1)*pt(1) + x(2)*pt(2) + x(3);
  
