function loadunload()

  ev_inputs = [0.0904 0.1808 0.2712 0.3616 0.452];
  ev_lambda = exp(-ev_inputs/3)

  ep_inputs = [0.06279244783307997 0.14509576209836975 0.2259340558658385 0.3019145714912858 0.3200070101047178]
  ep_lambda = exp(-ep_inputs/3)
  
  delT = 0.1;

  eps = 0.0;
  F0 = [1.000  0.0  0.0;   0.0  1.000  0.0;   0.0  0.0  1.000];
  F = [9.70316149436501e-01  0.0  0.0;   0.0  9.70316149436501e-01  0.0;   0.0  0.0  9.70316149436501e-01];
  [eps, lneps] = computeStrain(F0, F, eps, delT);
  lneps

  F0 = F;
  F = [9.79286713252325e-01  0.0  0.0;   0.0  9.79286713252325e-01  0.0;   0.0  0.0  9.79286713252325e-01];
  [eps, lneps] = computeStrain(F0, F, lneps, delT);
  lneps

  F0 = F;
  F = [9.41513429857278e-01  0.0  0.0;   0.0  9.41513429857278e-01  0.0;   0.0  0.0  9.41513429857278e-01];
  [eps, lneps] = computeStrain(F0, F, lneps, delT);
  lneps

  F0 = F;
  F = [9.52785714689278e-01  0.0  0.0;   0.0  9.52785714689278e-01  0.0;   0.0  0.0  9.52785714689278e-01];
  [eps, lneps] = computeStrain(F0, F, lneps, delT);
  lneps

  F0 = F;
  F = [9.13565685901867e-01  0.0  0.0;   0.0  9.13565685901867e-01  0.0;   0.0  0.0  9.13565685901867e-01];
  [eps, lneps] = computeStrain(F0, F, lneps, delT);
  lneps

  F0 = F;
  F = [9.27454676543165e-01  0.0  0.0;   0.0  9.27454676543165e-01  0.0;   0.0  0.0  9.27454676543165e-01];
  [eps, lneps] = computeStrain(F0, F, lneps, delT);
  lneps

  F0 = F;
  F = [8.86447538601615e-01  0.0  0.0;   0.0  8.86447538601615e-01  0.0;   0.0  0.0  8.86447538601615e-01];
  [eps, lneps] = computeStrain(F0, F, lneps, delT);
  lneps

  F0 = F;
  F = [9.04260143619469e-01  0.0  0.0;   0.0  9.04260143619469e-01  0.0;   0.0  0.0  9.04260143619469e-01];
  [eps, lneps] = computeStrain(F0, F, lneps, delT);
  lneps

  F0 = F;
  F = [8.60134362333383e-01  0.0  0.0;   0.0  8.60134362333383e-01  0.0;   0.0  0.0  8.60134362333383e-01];
  [eps, lneps] = computeStrain(F0, F, lneps, delT);
  lneps

  F0 = F;
  F = [8.98823131187731e-01  0.0  0.0;   0.0  8.98823131187731e-01  0.0;   0.0  0.0  8.98823131187731e-01];
  [eps, lneps] = computeStrain(F0, F, lneps, delT);
  lneps

  F0 = F;
  F = [1.001  0.0  0.0;   0.0  1.001  0.0;   0.0  0.0  1.001];
  [eps, lneps] = computeStrain(F0, F, lneps, delT)
  lneps
end

function [eps1, lneps1] = computeStrain(F0, F, eps, delT, ll)
  %Id = [1 0 0; 0 1 0; 0 0 1];
  %C = F'*F;
  %B = F*F';
  %E = trace(0.5*(C-Id))
  %e = trace(0.5*(Id-B))
  J = det(F);
  lnJ = log(J)

  nt = 10000;
  dt = delT/nt;
  dts = ones(nt)*dt;
  t = 0.0;
  lneps1 = 0;
  for i=1:length(dts)
    s = t/delT;
    Fn0 = (1-s)*F0 + s*F; 
    t = t + dt;
    s = t/delT;
    Fn1 = (1-s)*F0 + s*F; 
    Fdot = (Fn1 - Fn0)/dt;
    Finv = inv(Fn1);
    l = Fdot * Finv;
    d = 0.5 * (l + l');
    lnepsn1 = trace(d*dt);
    lneps1 = lneps1 + lnepsn1;
  end
  lneps1 = lneps1 + eps;
  eps1 = exp(lneps1)-1;
end
