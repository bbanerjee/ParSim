#include <gtest/gtest.h>
#include "../DruckerPragerInducedAniso.hpp"
#include <iostream>
#include <vector>

using namespace Vaango::Models::DruckerPragerInducedAniso;

extern "C" {
  void diamm_calc_(int* nblk, int* ninsv, double* dt, double* ui, 
                   double* stress, double* d, double* sv, double* usm);
}

TEST(DruckerPragerInducedAnisoTest, CompareWithFortran) {
    // 1. Material Properties
    MaterialProperties props;
    props.initialBulkModulus0 = 100.0e9;
    props.initialShearModulus0 = 75.0e9;
    props.yieldA1 = 200.0e6;
    props.yieldA2 = 0.1e-6;
    props.yieldA3 = 100.0e6;
    props.yieldA4 = 0.2;
    props.yieldA5 = 50.0e6;
    props.yieldA6 = 0.5;
    props.initialDensity = 1000.0;
    props.initialTemperature = 294.0;
    props.meltTemperature = 1000.0;
    props.homologousTempExponent = 1.0;
    props.useInducedAnisotropy = true;
    props.anisotropyCouplingG1 = 10.0e9;
    props.bulkModulusID = 1;
    props.bulkModulusParameter1 = 10.0e9;
    props.bulkModulusParameter2 = 1.0e9;

    // 2. Initial State
    StateVariables sv;
    sv.density = props.initialDensity;
    sv.temperature = props.initialTemperature;
    sv.jacobian = 1.0;
    sv.quasiStaticStress = {10.0e6, -5.0e6, 2.0e6, 1.0e6, 0.0, 0.0};
    sv.elasticStrain = {0.0001, -0.00005, 0.0, 0.0, 0.0, 0.0};

    // 3. Inputs
    double dt = 1.0e-6;
    SymmetricTensor strainInc = {0.0005, -0.0002, 0.0001, 0.0001, 0.0, 0.0};
    SymmetricTensor initialStress = sv.quasiStaticStress;

    // --- C++ Execution ---
    DruckerPragerInducedAnisoModel model(props);
    double usmCpp = 0.0;
    StateVariables svCpp = sv;
    model.updateStress(dt, strainInc, initialStress, svCpp, usmCpp);

    // --- Fortran Execution ---
    int nblk = 1;
    int ninsv = 26;
    std::vector<double> ui(47, 0.0);
    ui[0] = props.initialBulkModulus0;
    ui[1] = props.bulkModulusParameter1;
    ui[2] = props.bulkModulusParameter2;
    ui[3] = props.initialShearModulus0;
    ui[13] = props.useInducedAnisotropy ? 1.0 : 0.0;
    ui[4] = props.anisotropyCouplingG1 * 3.0; // G1_fortran = 3 * G1_cpp
    ui[7] = props.yieldA1 / (0.5773502691896257); // A1_fortran = A1_cpp / TOOR3
    ui[8] = props.yieldA2 / (1.7320508075688772); // A2_fortran = A2_cpp / ROOT3
    ui[9] = props.yieldA3 / (0.5773502691896257);
    ui[10] = props.yieldA4 * 1.7320508075688772; // A4_fortran = A4_cpp * ROOT3
    ui[14] = props.initialDensity;
    ui[15] = props.initialTemperature;
    ui[20] = props.meltTemperature;
    ui[25] = props.homologousTempExponent;
    ui[27] = (double)props.bulkModulusID;
    ui[28] = (double)props.shearModulusID;

    std::vector<double> svFort(ninsv, 0.0);
    svFort[4] = sv.volumetricStrain;
    svFort[5] = sv.temperature;
    svFort[7] = sv.density;
    svFort[9] = sv.jacobian;
    svFort[3] = sv.equivalentPlasticShearStrain;
    for(int i=0; i<6; ++i) {
        svFort[13+i] = sv.quasiStaticStress[i];
        svFort[19+i] = sv.elasticStrain[i];
    }

    double usmFort = 0.0;
    double stressFort[6];
    for(int i=0; i<6; ++i) stressFort[i] = initialStress[i];
    
    double dFort[6];
    for(int i=0; i<6; ++i) dFort[i] = strainInc[i] / dt;

    diamm_calc_(&nblk, &ninsv, &dt, ui.data(), stressFort, dFort, svFort.data(), &usmFort);

    // --- Comparison ---
    EXPECT_NEAR(usmCpp, usmFort, 1e-7);
    for(int i=0; i<6; ++i) {
        EXPECT_NEAR(svCpp.quasiStaticStress[i], svFort[13+i], 1e-2); // Stress scale is large
        EXPECT_NEAR(svCpp.elasticStrain[i], svFort[19+i], 1e-8);
    }
    EXPECT_NEAR(svCpp.equivalentPlasticShearStrain, svFort[3], 1e-8);
    EXPECT_NEAR(svCpp.temperature, svFort[5], 1e-6);
    EXPECT_NEAR(svCpp.density, svFort[7], 1e-6);
}
