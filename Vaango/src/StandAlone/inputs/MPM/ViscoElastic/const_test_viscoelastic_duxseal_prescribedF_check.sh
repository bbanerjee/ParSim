./partextract -partvar p.stress const_test_viscoelastic_duxseal_prescribedF.uda.000 > const_test_viscoelastic_duxseal_prescribedF.stress
#sig = load("const_test_viscoelastic_duxseal_prescribedF.stress");
#t = sig(:,1);
#sigx =  sig(:,5);
#plot(t, sigx);
