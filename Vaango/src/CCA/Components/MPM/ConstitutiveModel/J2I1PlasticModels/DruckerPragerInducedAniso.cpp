#include "DruckerPragerInducedAniso.hpp"
#include <iostream>

namespace Vaango::Models::DruckerPragerInducedAniso {

constexpr double ROOT3 = 1.7320508075688772;
constexpr double TOOR3 = 0.5773502691896257;
constexpr double TOOR2 = 0.7071067811865475;
constexpr double ROOT23 = 0.8164965809277260;

void DruckerPragerInducedAnisoModel::updateStress(double deltaTime,
                              const SymmetricTensor& strainInc,
                              const SymmetricTensor& initialStress,
                              StateVariables& sv,
                              double& uniaxialStrainModulus) 
{
    // 1. Initial Setup
    double totalStrainRateMag = getMagnitude(strainInc) / deltaTime;
    
    SymmetricTensor quasiStaticStressN;
    for(int i=0; i<6; ++i) {
        quasiStaticStressN[i] = sv.quasiStaticStress[i] * sv.jacobian;
    }

    double I1_n = getTrace(quasiStaticStressN);
    double rootJ2_n = getInvariantRootJ2(quasiStaticStressN);

    double tempMeltFactor = 1.0 - std::pow(std::max(sv.temperature - d_props.initialTemperature, 0.0) / 
                                          (d_props.meltTemperature - d_props.initialTemperature), 
                                          d_props.homologousTempExponent);

    // 2. Elastic Trial
    auto elasticStrainDevN = getDeviator(sv.elasticStrain);
    double beta2 = 0.5 * doubleDotProduct(elasticStrainDevN, elasticStrainDevN);
    double beta1 = getTrace(sv.elasticStrain);

    double K, G;
    calculateModuli(sv, beta1, beta2, sv.equivalentPlasticShearStrain, I1_n, rootJ2_n, K, G);
    
    uniaxialStrainModulus = (4.0 * G + 3.0 * K) / 3.0;
    double twoG = 2.0 * G;
    double twoG1 = 2.0 * d_props.anisotropyCouplingG1;

    double trDe = getTrace(strainInc);
    auto devDe = getDeviator(strainInc);

    SymmetricTensor trialStressInc;
    for(int i=0; i<6; ++i) {
        double delta = (i < 3) ? 1.0 : 0.0;
        trialStressInc[i] = K * trDe * delta + twoG * devDe[i];
        if (d_props.useInducedAnisotropy) {
             trialStressInc[i] += twoG1 * (trDe * elasticStrainDevN[i] + 
                                           doubleDotProduct(devDe, devDe) * delta);
        }
    }

    SymmetricTensor trialStress;
    for(int i=0; i<6; ++i) trialStress[i] = quasiStaticStressN[i] + trialStressInc[i];

    double lodeZT, lodeST;
    SymmetricTensor basisES;
    calculateLodeCoordinates(trialStress, lodeZT, lodeST, basisES);

    double yieldValue = calculateYieldFunction(lodeZT, lodeST, sv.equivalentPlasticShearStrain, tempMeltFactor);

    if (yieldValue <= 1e-2) {
        // --- Elastic Step ---
        for (int i = 0; i < 6; ++i) {
            sv.quasiStaticStress[i] = trialStress[i] / sv.jacobian;
            sv.elasticStrain[i] += strainInc[i];
        }
    } else {
        // --- Plastic Step (Newton Return Mapping) ---
        double z = lodeZT;
        double s = lodeST;

        for (int iter = 0; iter < 25; ++iter) {
            // Yield normal RN and flow direction RM
            double dFdz = (d_props.yieldA2 * d_props.yieldA3 * std::exp(d_props.yieldA2 * z) + d_props.yieldA4) * tempMeltFactor;
            double dFds = TOOR2;

            double dGdz = (d_props.yieldA2 * d_props.yieldA3 * std::exp(d_props.yieldA2 * z) + d_props.flowPotentialA4G) * tempMeltFactor;
            double dGds = TOOR2;

            // Flow direction components
            SymmetricTensor RM;
            for (int i = 0; i < 6; ++i) {
                double delta = (i < 3) ? 1.0 : 0.0;
                RM[i] = dGdz * TOOR3 * delta + dGds * basisES[i];
            }

            // Components of A tensor
            double RM_dot_EDEV = dGds * doubleDotProduct(basisES, elasticStrainDevN);
            double AZ = K * dGdz + ROOT3 * twoG1 * RM_dot_EDEV;
            double AS = twoG * dGds;

            // Scaled components of return direction tensor
            double fac = std::sqrt((z * z + s * s) / (dGdz * dGdz + dGds * dGds)) / K;
            double pz = fac * AZ;
            double ps = fac * AS;

            double beta = -yieldValue / (dFdz * pz + dFds * ps);

            z += beta * pz;
            s += beta * ps;

            if (std::abs(beta) < 1e-10) break;

            yieldValue = calculateYieldFunction(z, s, sv.equivalentPlasticShearStrain, tempMeltFactor);
        }

        // Updated quasistatic stress and increment
        SymmetricTensor qStressP;
        SymmetricTensor dTau;
        for (int i = 0; i < 6; ++i) {
            double delta = (i < 3) ? 1.0 : 0.0;
            qStressP[i] = z * TOOR3 * delta + s * basisES[i];
            dTau[i] = qStressP[i] - quasiStaticStressN[i];
        }

        double dz, ds;
        SymmetricTensor basisEDS;
        calculateLodeCoordinates(dTau, dz, ds, basisEDS);

        // Elastic strain increment
        double alpha1 = 0.0, alpha2 = 0.0;
        if (d_props.useInducedAnisotropy) {
            double zeta = (beta2 * twoG1 * twoG1) / G / K;
            double DDD = 6.0 * K * G * (1.0 - zeta) / (twoG1 * twoG1);
            double y1 = getTrace(dTau);
            double y2 = doubleDotProduct(elasticStrainDevN, dTau);
            alpha1 = (1.0 / DDD) * ((2.0 / 3.0) * beta2 * y1 / K - y2 / twoG1);
            alpha2 = (1.0 / DDD) * (1.5 * y2 / G - y1 / twoG1);
        }

        SymmetricTensor dEe;
        for (int i = 0; i < 6; ++i) {
            double delta = (i < 3) ? 1.0 : 0.0;
            dEe[i] = (1.0 / (3.0 * K)) * dz * TOOR3 * delta + (1.0 / twoG) * ds * basisEDS[i] +
                     alpha1 * delta + alpha2 * elasticStrainDevN[i];
            sv.elasticStrain[i] += dEe[i];
        }

        // Plastic strain increment
        SymmetricTensor dEp;
        for (int i = 0; i < 6; ++i) {
            dEp[i] = strainInc[i] - dEe[i];
        }
        auto dEpDev = getDeviator(dEp);
        double deqps = std::sqrt(doubleDotProduct(dEpDev, dEpDev) + (1.0 / 3.0) * std::pow(getTrace(dEp), 2));
        if (deqps < 1e-16) {
            deqps = 0.0;
        }
        sv.equivalentPlasticShearStrain += deqps;
        sv.plasticVolumeStrain += getTrace(dEp);

        for (int i = 0; i < 6; ++i) {
            sv.quasiStaticStress[i] = qStressP[i] / sv.jacobian;
        }
    }

    // 3. Final Updates (Temperature, Density, etc.)
    sv.volumetricStrain += trDe;
    sv.density *= std::exp(-trDe);
    sv.jacobian = d_props.initialDensity / sv.density;

    sv.stressInvariantI1 = getTrace(sv.quasiStaticStress);
    sv.stressInvariantRootJ2 = getInvariantRootJ2(sv.quasiStaticStress);
    sv.totalStrainRateMagnitude = totalStrainRateMag;
    sv.soundSpeed = std::sqrt(K / d_props.initialDensity);
    
    auto finalDevE = getDeviator(sv.elasticStrain);
    sv.anisotropyMeasure = calculateAnisotropyMeasure(K, G, d_props.anisotropyCouplingG1, 
                                                      doubleDotProduct(finalDevE, finalDevE));
}

void DruckerPragerInducedAnisoModel::calculateModuli(const StateVariables& sv, double ev, double beta2, 
                                 double eqps, double ri1, double rtj2,
                                 double& bulkModulus, double& shearModulus) {
    // Bulk Modulus Logic (IDK == 1)
    if (d_props.bulkModulusID == 1) {
        double p = -ri1 / 3.0;
        bulkModulus = d_props.initialBulkModulus0 + d_props.bulkModulusParameter1 * std::exp(-d_props.bulkModulusParameter2 / std::max(p, 1e-40));
    } else {
        bulkModulus = d_props.initialBulkModulus0;
    }

    // Shear Modulus Logic (IDG == 0 - SGC formula)
    if (d_props.shearModulusID == 0) {
        double rjFactor = std::pow(d_props.initialDensity / sv.density, 1.0/3.0);
        shearModulus = std::max(d_props.initialShearModulus0, 
                                d_props.initialShearModulus0 + d_props.shearModulusParameter1 * ev * rjFactor + 
                                d_props.shearModulusParameter2 * (sv.temperature - 300.0));
    } else {
        shearModulus = d_props.initialShearModulus0;
    }
}

double DruckerPragerInducedAnisoModel::calculateYieldFunction(double lodeZ, double lodeS, double eqps, double tempMeltFactor) const {
    double pressure = -lodeZ * TOOR3;
    double I1 = -3.0 * pressure;
    double flowStress = d_props.yieldA1 - d_props.yieldA3 * std::exp(d_props.yieldA2 * I1) - 
                        d_props.yieldA4 * I1 + d_props.yieldA5 * std::pow(ROOT23 * eqps, d_props.yieldA6);
    return TOOR2 * lodeS - flowStress * tempMeltFactor;
}

void DruckerPragerInducedAnisoModel::calculateLodeCoordinates(const SymmetricTensor& stress, double& z, double& s, SymmetricTensor& basisES) const {
    z = TOOR3 * getTrace(stress);
    auto dev = getDeviator(stress);
    s = getMagnitude(dev);
    double mag = (s < 1e-10) ? 1e90 : s;
    for(int i=0; i<6; ++i) basisES[i] = dev[i] / mag;
}

double DruckerPragerInducedAnisoModel::calculateAnisotropyMeasure(double K, double G, double G1, double devE2) const {
    double isoMag = std::sqrt(9.0 * K * K + 5.0 * (4.0 * G * G));
    double totalMag = std::sqrt(isoMag * isoMag + 6.0 * (4.0 * G1 * G1) * devE2);
    return (2.0 / M_PI) * std::acos(isoMag / totalMag);
}

} // namespace Vaango::Models::DruckerPragerInducedAniso
