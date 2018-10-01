function computeDefGrad()

  filename = "TabularCapTest_11_PrescribedDeformation.inp";
  fid = fopen(filename, "w");
  t0 = 0;
  t1 = 0.75; 
  t2 = 1.0;
  E0 = [0 0 0; 0 0 0; 0 0 0];
  E1 = [1.5 0 0; 0 -0.75 0; 0 0 -0.75];
  E2 = [1 0 0; 0 0.25 0; 0 0 0.25];

  dt = 0.0001;

  F0 = expm(E0);
  fprintf(fid, "%f %.16f %f %f %f %.16f %f %f %f %.16f 0 0 0 1\n", ...
    t0, F0(1,1), F0(1,2), F0(1,3), ...
    F0(2,1), F0(2,2), F0(2,3), ...
    F0(3,1), F0(3,2), F0(3,3));
  F01 = computeF(F0, E0, E1, t0, t1, dt, fid);
  F1 = F01{length(F01), 1};
  F12 = computeF(F1, E1, E2, t1, t2, dt, fid);
  F12
  fclose(fid);
end

function [F] = computeF(F0, E0, E1, t0, t1, dt, fid)

  Edot = (E1 - E0)/(t1 - t0);

  nt = ceil((t1 - t0)/dt);
  times = linspace(t0, t1, nt);

  F = cell(nt, 1);
  F(1,1) = F0;
  for ii = 2:length(times)
    Ldt = Edot * (times(ii) - t0);
    expLdt = expm(Ldt);
    F1 = expLdt * F0;
    F(ii, 1) = F1;
    fprintf(fid, "%f %.16f %f %f %f %.16f %f %f %f %.16f 0 0 0 1\n", ...
      times(ii), F1(1,1), F1(1,2), F1(1,3), ...
      F1(2,1), F1(2,2), F1(2,3), ...
      F1(3,1), F1(3,2), F1(3,3));
  end
end
