#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_BorjaT.h>
#include <iostream>

#include <gtest/gtest.h>

TEST(ModelStateCRTPTest, Constructors)
{
  Vaango::ModelState_BorjaT state;
  state.porosity = 1.0;
  state.pressure = 3.14e6;
  //std::cout << "p = " << state.pressure << " phi = " << state.porosity << "\n";

  auto state_copy(state);
  //std::cout << "p copy = " << state_copy.pressure << " phi copy = " << state_copy.porosity << "\n";
  ASSERT_DOUBLE_EQ(state.pressure, state_copy.pressure);

  auto state_pcopy(&state);
  //std::cout << "p pcopy = " << state_pcopy->pressure << " phi pcopy = " << state_pcopy->porosity << "\n";
  ASSERT_DOUBLE_EQ(state.pressure, state_pcopy->pressure);

  auto state_ecopy = state;
  //std::cout << "p ecopy = " << state_ecopy.pressure << " phi ecopy = " << state_ecopy.porosity << "\n";
  ASSERT_DOUBLE_EQ(state.pressure, state_ecopy.pressure);

  auto state_epcopy = &state;
  //std::cout << "p epcopy = " << state_epcopy->pressure << " phi epcopy = " << state_epcopy->porosity << "\n";
  ASSERT_DOUBLE_EQ(state.pressure, state_epcopy->pressure);

  Uintah::Matrix3 stress(100, 200, 300, 200, 400, 500, 300, 500, 600);
  auto s = state_ecopy.updateStressInvariants(stress);
  //std::cout << std::setprecision(16)
  //          << "I1 orig = " << state.I1 << " J2 orig = " << state.J2 << " p orig = " << state.p 
  //          << " q orig = " << state.q << "\n";
  //std::cout << std::setprecision(16)
  //          << "I1 ecopy = " << state_ecopy.I1 << " J2 ecopy = " << state_ecopy.J2 
  //          << " p ecopy = " << state_ecopy.p << " q ecopy = " << state_ecopy.q << "\n";
  //std::cout << "dev(stress) = " << s << "\n";
  ASSERT_NEAR(s.Trace(), 0.0, 1.0e-12);
  ASSERT_DOUBLE_EQ(state_ecopy.I1, 1100.0);
  ASSERT_DOUBLE_EQ(state_ecopy.J2, 443333.33333333337);
  ASSERT_NEAR(state_ecopy.p, 366.6666666666667, 1.0e-8);
  ASSERT_NEAR(state_ecopy.q, 1153.2562594670796, 1.0e-8);

  state_ecopy.updateStressInvariants();
  //std::cout << std::setprecision(16)
  //          << "I1 ecopy = " << state_ecopy.I1 << " J2 ecopy = " << state_ecopy.J2 
  //          << " p ecopy = " << state_ecopy.p << " q ecopy = " << state_ecopy.q << "\n";
  ASSERT_NEAR(state_ecopy.q, 1153.2562594670796, 1.0e-8);

  auto num = state_ecopy.numStateVar();
  //std::cout << "num = " << num << "\n";
  ASSERT_EQ(num, 33);

  auto strain = Uintah::Matrix3(1, 2, 3, 2, 4, 5, 3, 5, 6);
  auto strain_tr = Uintah::Matrix3(2, 3, 4, 3, 5, 6, 4, 6, 7);
  auto dev_strain = state.updateStrainScalars(strain, strain_tr);
  state.p_c = 1000.0;
  state.p_c0 = 800.0;
  //std::cout << std::setprecision(16)
  //          << "eps_dev = " << dev_strain.first << "\n" << dev_strain.second << "\n";
  //std::cout << " epse_v = " << state.epse_v << " epse_v_tr = " << state.epse_v_tr 
  //          << " epse_s = " << state.epse_s << " epse_s_tr = " << state.epse_s_tr << "\n";
  //std::cout << " p_c = " << state_epcopy->p_c << " p_c0 = " << state_epcopy->p_c0 << "\n";
  ASSERT_EQ(state.epse_s, 7.688375063113864);
  ASSERT_EQ(state.epse_s_tr, 9.47511360236793);
}

