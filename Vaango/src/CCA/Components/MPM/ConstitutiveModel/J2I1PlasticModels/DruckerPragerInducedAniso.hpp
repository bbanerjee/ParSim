#ifndef VAANGO_CCA_COMPONENTS_MPM_CONSTITUTIVEMODEL_J2I1PLASTICMODELS_DIAMM_HPP
#define VAANGO_CCA_COMPONENTS_MPM_CONSTITUTIVEMODEL_J2I1PLASTICMODELS_DIAMM_HPP

#include <array>
#include <cmath>
#include <vector>
#include <optional>
#include <algorithm>
#include <span>

namespace Vaango::Models::DruckerPragerInducedAniso {

// Standard 2nd order symmetric tensor (ordered: 11, 22, 33, 12, 23, 13)
using SymmetricTensor = std::array<double, 6>;

struct MaterialProperties {
    double initialBulkModulus0{0.0};
    double bulkModulusParameter1{0.0};
    double bulkModulusParameter2{0.0};
    double initialShearModulus0{0.0};
    double shearModulusParameter1{0.0};
    double shearModulusParameter2{0.0};
    double shearModulusParameter3{0.0};
    double anisotropyCouplingG1{0.0};
    
    double yieldA1{0.0};
    double yieldA2{0.0};
    double yieldA3{0.0};
    double yieldA4{0.0};
    double yieldA5{0.0};
    double yieldA6{0.0};
    double flowPotentialA4G{0.0};

    double initialDensity{0.0};
    double initialTemperature{0.0};
    double initialSoundSpeed{0.0};
    double gruneisenParameter{0.0};
    double specificHeat{0.0};
    double meltTemperature{0.0};
    double homologousTempExponent{0.0};
    double taylorQuinneyCoeff{1.0};

    bool useInducedAnisotropy{false};
    bool isRateDependent{false};
    
    // Rate dependence (Overstress) parameters
    double tau1{0.0}, tau2{0.0}, tau3{0.0}, tau4{0.0};

    // EOS IDs
    int bulkModulusID{0};
    int shearModulusID{0};
};

struct StateVariables {
    double totalStrainRateMagnitude{0.0};
    double stressInvariantI1{0.0};
    double stressInvariantRootJ2{0.0};
    double equivalentPlasticShearStrain{0.0};
    double volumetricStrain{0.0};
    double temperature{0.0};
    double soundSpeed{0.0};
    double density{0.0};
    double internalEnergy{0.0};
    double jacobian{1.0};
    double anisotropyMeasure{0.0};
    double plasticVolumeStrain{0.0};
    
    SymmetricTensor quasiStaticStress{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    SymmetricTensor elasticStrain{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
};

// --- Math Utilities ---

inline double getTrace(const SymmetricTensor& A) {
    return A[0] + A[1] + A[2];
}

inline double doubleDotProduct(const SymmetricTensor& A, const SymmetricTensor& B) {
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + 
           2.0 * (A[3]*B[3] + A[4]*B[4] + A[5]*B[5]);
}

inline double getMagnitude(const SymmetricTensor& A) {
    return std::sqrt(doubleDotProduct(A, A));
}

inline SymmetricTensor getDeviator(const SymmetricTensor& A) {
    double mean = getTrace(A) / 3.0;
    return {A[0] - mean, A[1] - mean, A[2] - mean, A[3], A[4], A[5]};
}

inline double getInvariantRootJ2(const SymmetricTensor& A) {
    auto dev = getDeviator(A);
    return std::sqrt(0.5 * doubleDotProduct(dev, dev));
}

// --- Model Logic ---

class DruckerPragerInducedAnisoModel {
public:
    explicit DruckerPragerInducedAnisoModel(const MaterialProperties& props) : d_props(props) {}

    void updateStress(double deltaTime,
                      const SymmetricTensor& strainIncrement,
                      const SymmetricTensor& initialStress,
                      StateVariables& sv,
                      double& uniaxialStrainModulus);

private:
    MaterialProperties d_props;

    void calculateModuli(const StateVariables& sv, 
                         double elasticVolStrain, 
                         double elasticDevInv2,
                         double eqps,
                         double I1, double rootJ2,
                         double& bulkModulus, 
                         double& shearModulus);

    double calculateYieldFunction(double lodeZ, double lodeS, 
                                  double eqps, double tempMeltFactor) const;

    void calculateLodeCoordinates(const SymmetricTensor& stress, 
                                  double& z, double& s, 
                                  SymmetricTensor& basisES) const;

    double calculateAnisotropyMeasure(double K, double G, double G1, double devE2) const;
};

} // namespace Vaango::Models::DruckerPragerInducedAniso

#endif
