function testMGE()

  J = linspace(0.001, 1.1, 100);
  p = [];
  dp_dJ = [];
  for i = 1:length(J)
    p(i) = pressure(J(i));
    dp_dJ(i) = deriv_p_J(J(i));
  end

  figure;
  plot(J, p, '-', 'linewidth', 2); hold on;
  grid on;

  figure;
  plot(J, dp_dJ, '-', 'linewidth', 2); hold on;
  grid on;

  eta = 1. - J;
  eta2 = eta .* eta;
  eta3 = eta2 .* eta;
  S_1 = 1.123;
  S_2 = 3.98;
  S_3 = -5.8;
  alpha = (1. - S_1 * eta - S_2 * eta2 - S_3 * eta3);
  alpha2 = alpha .* alpha;

  J_min_plus = 1 + (S_2 + sqrt(S_2 * S_2 - 3 * S_1 * S_3))/ (3 * S_3);
  J_min_minus = 1 + (S_2 - sqrt(S_2 * S_2 - 3 * S_1 * S_3))/ (3 * S_3);
  figure;
  plot(J, alpha2, '-', 'linewidth', 2); hold on;
  plot(J_min_plus, 0.0, 'x', 'linewidth', 2);
  plot(J_min_minus, 0.0, 'o', 'linewidth', 2);
  grid on;

  dp_dJ_0 = deriv_p_J(1.0)

end

function p = pressure(J)

  rho_0 = 2200;
  C_0 = 1680;
  rho0_C0_sq = rho_0 * C_0 * C_0;
  S_1 = 1.123;
  S_2 = 3.98;
  S_3 = -5.8;
  Gamma_0 = 0.59;
  e = 0.0;

  eta = 1. - J;

  p = rho_0 * Gamma_0 * e;
  if (eta >= 0.0) 
    eta2 = eta * eta;
    eta3 = eta2 * eta;
    alpha = (1. - S_1 * eta - S_2 * eta2 - S_3 * eta3);
    alpha2 = alpha * alpha; 
    p = p + rho0_C0_sq * eta * (1. - .5 * Gamma_0 * eta) / alpha2;
  else 
    p = p + rho0_C0_sq * eta;
  end

end

function dp_dJ = deriv_p_J(J)

  rho_0 = 2200;
  C_0 = 1680;
  rho0_C0_sq = rho_0 * C_0 * C_0;
  S_1 = 1.123;
  S_2 = 3.98;
  S_3 = -5.8;
  Gamma_0 = 0.59;

  eta = 1. - J;

  dp_dJ = 0.0;
  if (eta >= 0.0) 
    eta2 = eta * eta;
    eta3 = eta2 * eta;
    alpha = 1. - S_1 * eta - S_2 * eta2 - S_3 * eta3;
    dalpha_deta = - S_1 - 2.0 * S_2 * eta - 3.0 * S_3 * eta2;
    alpha2 = alpha * alpha;
    dp_dJ = - (rho0_C0_sq / alpha2) * ...
              ( 1 - Gamma_0 * eta  - ...
                2.0 * eta * (1 - Gamma_0/2 * eta) / alpha * dalpha_deta);
  else 
    dp_dJ = - rho0_C0_sq;
  end

end
