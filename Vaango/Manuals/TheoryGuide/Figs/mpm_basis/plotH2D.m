function plotH2D()

  xp = 0.0;
  yp = 0.0;
  lp = 0.25;

  numPts = 200;
  y = linspace(-1, 1, numPts);
  x = linspace(-1, 1, numPts);

  Pxy = Pulse2D(x, y, xp, yp, lp);
  %plot3(x, y, Pxy);

  xg = 0.0;
  yg = 0.0;
  lg = 0.25;
  Sxy = Linear2DXY(x, y, xg, yg, lg);
  %plot3(x, y, Sxy);
  mesh(x, y, Sxy);
  
end
  
function Pxy = Pulse2D(x, y, xp, yp, lp)
  Px = Pulse(x, xp, lp);
  Py = Pulse(y, yp, lp);
  Pxy = Px' .* Py;
end

function Px = Pulse(x, xp, lp)
  Px_minus = Heaviside(x, xp-lp);
  Px_plus = Heaviside(x, xp+lp);
  Px = Px_minus - Px_plus;
end
  
function Hxy = Heaviside2D(x, y, xp, yp)
  Hx = Heaviside(x, xp);
  Hy = Heaviside(y, yp);
  Hxy = Hx' .* Hy;
end

function Hx = Heaviside(x, xp)
  Hx = zeros(size(x));
  Hx(x >= xp) = 1;
end

function Sxy = Linear2DXY(x, y, xg, yg, lg)
  Sx = LinearX(x, xg, lg);
  Sy = LinearX(y, yg, lg);
  Sxy = Sx' .* Sy;
end

function Sx = LinearX(x, xg, lg)
  xi = zeros(size(x));
  xig = zeros(size(x));
  xi(x < xg) = Xi(x(x < xg), xg-lg, xg);
  xi(x >= xg) = Xi(x(x >= xg), xg, xg+lg);
  xig(x < xg) = -1;
  xig(x >= xg) = 1;
  Sx = Linear(xi, xig);
end

function xi = Xi(x, x1, x2)
  xi = 2.0 * (x - x1)/(x2 - x1) - 1;
end

function Sxi = Linear(xi, xig)
  Sxi = 0.5 * (1 - xi .* xig);
  Sxi(xi < -1.0 | xi > 1.0) = 0;
end
