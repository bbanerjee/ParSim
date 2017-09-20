minEps = -10;
maxEps = 10;
Eps = linspace(minEps, maxEps, 30);
Sig = 3000 .* tanh(0.25*sigmoid_x)
plot(Eps, Sig, 'k-', 'LineWidth', 2); hold on;

PlasticStrainVol =  [-2.5, 2.5];
TotalStrainVol(1,:) =  Eps + PlasticStrainVol(1);
Pressure(1,:) = Sig;
TotalStrainVol(2,:) =  Eps + PlasticStrainVol(2);
Pressure(2,:) = Sig;

ElasticStrainVol = TotalStrainVol - PlasticStrainVol';
plot(TotalStrainVol(1,:), Pressure(1,:), 'b-'); hold on;
plot(TotalStrainVol(2,:), Pressure(2,:), 'b-'); hold on;
plot(ElasticStrainVol(1,:), Pressure(1,:), 'r-'); hold on;
plot(ElasticStrainVol(2,:), Pressure(2,:), 'r-'); hold on;
PlasticStrainVol =  [-5.0, 5.0];
TotalStrainVol = ElasticStrainVol + PlasticStrainVol';
plot(TotalStrainVol(1,:), Pressure(1,:), 'g-', 'LineWidth', 2); hold on;
plot(TotalStrainVol(2,:), Pressure(2,:), 'g-', 'LineWidth', 2); hold on;
printf('%f, ', TotalStrainVol(1,:))
printf('\n')
printf('%f, ', TotalStrainVol(2,:))
printf('\n')
printf('%f, ', Pressure(2,:))
printf('\n')
