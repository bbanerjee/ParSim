function [theta_bisect, theta_newton] = ellipticEinv(phi_max, s_tilde, mm, theta_min, theta_max, tol)

  theta_bisect = ellipticEinvBisect(phi_max, s_tilde, mm, theta_min, theta_max, tol);
  %% theta_newton = ellipticEinvNewton(s_tilde, mm, tol);
  theta_newton = theta_bisect;

%
% Find inverse of elliptic integral E for fixed m
% Bisection method
function [theta] = ellipticEinvBisect(phi_max, s_tilde, mm, theta_min, theta_max, tol)

  if nargin<4, tol = eps; end
  if nargin<3, error('Not enough input arguments.'); end

  if ~isreal(s_tilde) || ~isreal(mm)
    error('Input arguments must be real.');
  end

  if length(mm)==1, mm = mm(ones(size(s_tilde))); end
  if length(s_tilde)==1, s_tilde = s_tilde(ones(size(mm))); end
  if ~isequal(size(mm),size(s_tilde)), error('s_tilde and M must be the same size.'); end

  theta = zeros(1,length(s_tilde)); 
  for ii=1:length(s_tilde)

    % Save m and E
    m_in = mm(ii);
    s_tilde_in = s_tilde(ii);

    % Special cases
    if (abs(s_tilde_in) < eps)
      theta(ii) = 0.0;
      continue;
    end
    if (m_in == 1.0 && s_tilde_in == 1.0)
      theta(ii) = pi/2.0;
      continue;
    end

    count_max = 100;
    count = 1;
    theta_mid_k = -1.0;
    theta_min_k = theta_min*pi/phi_max;
    theta_max_k = theta_max*pi/phi_max;
    while (count < count_max)
      theta_mid_k = 0.5*(theta_min_k + theta_max_k);
      [ff_mid, ~] = computeF(theta_mid_k, m_in, s_tilde_in);
      if (abs(ff_mid) < tol || (theta_max_k-theta_min_k)*0.5 < tol)
         break;
      end
      count = count + 1;
      [ff_min, ~] = computeF(theta_min_k, m_in, s_tilde_in);
      if (sign(ff_mid) == sign(ff_min))
        theta_min_k = theta_mid_k;
      else
        theta_max_k = theta_mid_k;
      end
    end
    if (count == count_max)
       ['Did not converge']
    end
    theta(ii) = theta_mid_k;
  end

% Newton method
function [theta] = ellipticEinvNewton(s_tilde, mm, tol)

  if nargin<4, tol = eps; end
  if nargin<3, error('Not enough input arguments.'); end

  if ~isreal(s_tilde) || ~isreal(mm)
    error('Input arguments must be real.');
  end

  if length(mm)==1, mm = mm(ones(size(s_tilde))); end
  if length(s_tilde)==1, s_tilde = s_tilde(ones(size(mm))); end
  if ~isequal(size(mm),size(s_tilde)), error('s_tilde and M must be the same size.'); end

  theta = zeros(1,length(s_tilde)); 
  for ii=1:length(s_tilde)

    % Save m and E
    m_in = mm(ii);
    s_tilde_in = s_tilde(ii);

    % Special cases
    if (abs(s_tilde_in) < eps)
      theta(ii) = 0.0;
      continue;
    end
      
    if (m_in == 1.0 && s_tilde_in == 1.0)
      theta(ii) = pi/2.0;
      continue;
    end

    % inputs
    zz = s_tilde_in; m_c = 1.0 - m_in;

    % complete integral initialization
    [~,E1] = ellipke(m_in, tol); 

    zeta = 1.0 - zz./E1;
    rr = sqrt(zeta.*zeta+m_c.*m_c); 
    theta_m = atan(m_c./(zz+eps));

    % “Empirical” initialization [1]
    theta_k = pi/2.0 + sqrt(rr).*(theta_m - (pi/2)); 
    [ff, dfdtheta] = computeF(theta_k, m_in, s_tilde_in);

    count = 0;
    while (abs(ff) > tol)
      theta_k = theta_k - ff/dfdtheta;
      [ff, dfdtheta] = computeF(theta_k, m_in, s_tilde_in);
      count = count + 1;
      if (count > 10)
        [theta_k m_in s_tilde_in ff dfdtheta tol]
        'Did not converge'
        break;
      end
    end

    % Compute theta
    theta(ii) = theta_k;
  end

% Compute the residual (f)
function [ff, dfdtheta] = computeF(theta, mu, stilde)

  % Compute incomplete elliptic integral
  [~,EE_i] = elliptic123(theta, mu);

  % Compute f
  ff = stilde - EE_i;

  % Compute dfdtheta
  dfdtheta = -sqrt(1.0 - mu*sin(theta)^2);


