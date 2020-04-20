#include "StateMohrCoulomb.h"

#include <iostream>
#include <cmath>

using namespace Uintah;

StateMohrCoulomb::StateMohrCoulomb()
{
  stress = Vector6::Zero(); 
  strain = Vector7::Zero(); 
  plasticStrain = Vector6::Zero(); 
  microStress = Vector6::Zero();
  microStrain = Vector7::Zero();
  microPlasticStrain = Vector6::Zero();
  state = Vector3::Zero();
  microState = Vector3::Zero();
  state(0) = 1e6; // p0Star, unused
  state(1) = 1e8; // Yield suction, unused
  state(2) = 1.5; // Specific Volume, unused
}

void
StateMohrCoulomb::update(const Vector6& plasticStrainInc, 
                         const Vector7& strainInc,
                         const Vector6& stressInc, 
                         double p0StarInc)
{
  double microcheck;
  for (int i = 0; i < 7; i++) {
    microcheck = strain(i);
    microStrain(i) += strainInc(i);
    strain(i) += microStrain(i);
    microStrain(i) -= (strain(i) - microcheck);
  }

  for (int i = 0; i < 6; i++) {
    microcheck = stress(i);
    microStress(i) += stressInc(i);
    stress(i) += microStress(i);
    microStress(i) -= (stress(i) - microcheck);
  }

  for (int i = 0; i < 6; i++) {
    microcheck = plasticStrain(i);
    microPlasticStrain(i) += plasticStrainInc(i);
    plasticStrain(i) += microPlasticStrain(i);
    microPlasticStrain(i) -= (plasticStrain(i) - microcheck);
  }

  microcheck = state(0);
  microState(0) += p0StarInc;
  state(0) += microState(0);
  microState(0) -= (state(0) - microcheck);

  microcheck = state(2);
  microState(2) += state(2) * (exp(-(strainInc(0) + strainInc(1) + strainInc(2))) - 1);
  // microState(2) -= (strainInc(0)+strainInc(1)+strainInc(2))*state(2);
  state(2) += microState(2); // updated specific volume
  microState(2) -= (state(2) - microcheck);
}

double
StateMohrCoulomb::meanStress() const
{
  return firstInvariant() / 3.0;
}

double
StateMohrCoulomb::shearStress() const
{
  double J2 = secondDevInvariant();
  double shearStress = 0.0;
  if (J2 > 1.0-14) {
    shearStress = std::sqrt(3.0 * J2); 
  }
  return shearStress;
}

double
StateMohrCoulomb::firstInvariant() const
{
  double I1 =  (stress(0) + stress(1) + stress(2));
  return I1;
}

double
StateMohrCoulomb::secondInvariant() const
{
  double I2 = stress(0) * stress(1) + stress(1) * stress(2) +
              stress(2) * stress(0) - stress(3) * stress(3) -
              stress(4) * stress(4) - stress(5) * stress(5);
  return I2;
}

double
StateMohrCoulomb::thirdInvariant() const
{
  double I3 =
    stress(0) * stress(1) * stress(2) + 2 * stress(3) * stress(4) * stress(5) -
    stress(0) * stress(5) * stress(5) -     stress(1) * stress(4) * stress(4) -
    stress(2) * stress(3) * stress(3);
  return I3;
}

double
StateMohrCoulomb::firstDevInvariant() const
{
  return 0.0;
}

double
StateMohrCoulomb::secondDevInvariant() const
{
  double J2 = ((stress(0) - stress(1)) * (stress(0) - stress(1)) +
               (stress(0) - stress(2)) * (stress(0) - stress(2)) +
               (stress(1) - stress(2)) * (stress(1) - stress(2))) / 6.0 +
               (stress(3) * stress(3) + stress(4) * stress(4) + 
                stress(5) * stress(5));
  return J2;
}

double
StateMohrCoulomb::thirdDevInvariant() const
{
  double I1 = firstInvariant();
  double I2 = secondInvariant();
  double I3 = thirdInvariant();
  double J3 = I1 * I1 * I1 * 2.0 / 27.0 - I1 * I2 / 3.0 + I3;

  return J3;
}

double
StateMohrCoulomb::getTheta()
{
  double J2 = secondDevInvariant();
  double J3 = thirdDevInvariant();
  double sin3Theta = 0.0;
  if (J2 != 0) {
    sin3Theta = std::sqrt(3.0) * (-1.5) * J3 / std::pow(J2, 1.5);
  }
  if (sin3Theta > 1) {
    sin3Theta = 1.0;
  } else if (sin3Theta < -1) {
    sin3Theta = -1.0;
  }
  double theta = (std::asin(sin3Theta)) / 3.0;

  return theta;
}

double
StateMohrCoulomb::getThetaDeg()
{
  double thetaDeg = getTheta();
  thetaDeg *= (45.0 / std::atan(1.0));
  if (thetaDeg < 0) {
    thetaDeg += 360.0;
  } else if (thetaDeg > 360) {
    thetaDeg -= 360.0;
  }

  return thetaDeg;
}

double
StateMohrCoulomb::getThetaDeg_0()
{
  double thetaDeg = getTheta();
  thetaDeg *= (45.0 / std::atan(1.0));
  thetaDeg -= 30;
  if (thetaDeg < 0) {
    thetaDeg += 360.0;
  } else if (thetaDeg > 360) {
    thetaDeg -= 360.0;
  }

  return thetaDeg;
}

bool
StateMohrCoulomb::checkIfFinite() const
{
  bool F = true;
  for (int i = 0; i < 6; i++)
    if (!std::isfinite(stress(i))) {
      F = false;
      std::cout << "Stress[" << i << "]=" << stress(i) << "\n";
      std::cout << "_finite(Stress[" << i << "])=" << std::isfinite(stress(i)) << "\n";
    }
  for (int i = 0; i < 7; i++)
    if (!std::isfinite(strain(i))) {
      F = false;
      std::cout << "Strain[" << i << "]=" << strain(i) << "\n";
      std::cout << "_finite(Strain[" << i << "])=" << std::isfinite(strain(i)) << "\n";
    }
  for (int i = 0; i < 3; i++)
    if (!std::isfinite(state(i))) {
      F = false;
      std::cout << "State[" << i << "]=" << state(i) << "\n";
      std::cout << "_finite(State[" << i << "])=" << std::isfinite(state(i)) << "\n";
    }
  if (!F) {
    std::cout << "Point internal values incorrect." << "\n";
    std::cout << "Stress:" << stress(0) << " " << stress(1) << " " << stress(2)
         << " " << stress(3) << " " << stress(4) << " " << stress(5) << " "
         << "\n";
    std::cout << "Strain:" << strain(0) << " " << strain(1) << " " << strain(2)
         << " " << strain(3) << " " << strain(4) << " " << strain(5) << " "
         << strain(6) << "\n";
    std::cout << "State variables:" << state(0) << " " << state(1) << " " << state(2)
         << "\n";
    std::cout << "Stress finite:" << std::isfinite(stress(0)) << " "
         << std::isfinite(stress(1)) << " " << std::isfinite(stress(2)) << " "
         << std::isfinite(stress(3)) << " " << std::isfinite(stress(4)) << " "
         << std::isfinite(stress(5)) << " " << "\n";
    std::cout << "strain finite:" << std::isfinite(strain(0)) << " "
         << std::isfinite(strain(1)) << " " << std::isfinite(strain(2)) << " "
         << std::isfinite(strain(3)) << " " << std::isfinite(strain(4)) << " "
         << std::isfinite(strain(5)) << " " << "\n";
    std::cout << "State variables finite:" << std::isfinite(state(0)) << " "
         << std::isfinite(state(1)) << " " << std::isfinite(state(2)) << "\n";
  }

  return F;
}

void 
StateMohrCoulomb::setStressEigen(const Vector3& eigenvals)
{
  stress(0) = eigenvals[2];
  stress(1) = eigenvals[1];
  stress(2) = eigenvals[0];
  stress(3) = 0.0;
  stress(4) = 0.0;
  stress(5) = 0.0;

}
