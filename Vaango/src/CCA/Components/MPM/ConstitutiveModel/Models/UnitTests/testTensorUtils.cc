#include <CCA/Components/MPM/ConstitutiveModel/Utilities/TensorUtils.h>

#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include <gtest/gtest.h>

using namespace Vaango;

TEST(TensorUtilsTest, constructors)
{
  Uintah::Matrix3 I(1, 0, 0, 0, 1, 0, 0, 0, 1);
  Tensor::Vector6Mandel I_vec = Tensor::constructVector6Mandel(I);  
  Tensor::Vector6Mandel I_vec_test; 
  I_vec_test << 1, 1, 1, 0, 0, 0;
  for (auto ii = 0u; ii < 6; ++ii) {
    ASSERT_DOUBLE_EQ(I_vec(ii), I_vec_test(ii));
  }

  Tensor::Matrix6Mandel I_mat = Tensor::constructMatrix6Mandel(I, I);  
  Tensor::Matrix6Mandel I_mat_test; 
  I_mat_test << 1, 1, 1, 0, 0, 0,
                1, 1, 1, 0, 0, 0,
                1, 1, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0;
  for (auto ii = 0u; ii < 6; ++ii) {
    for (auto jj = 0u; jj < 6; ++jj) {
      ASSERT_DOUBLE_EQ(I_mat(ii, jj), I_mat_test(ii, jj));
    }
  }

  Tensor::Matrix6Mandel I_mat_alt = Tensor::constructMatrix6Mandel(I_vec, I_vec);  
  for (auto ii = 0u; ii < 6; ++ii) {
    for (auto jj = 0u; jj < 6; ++jj) {
      ASSERT_DOUBLE_EQ(I_mat_alt(ii, jj), I_mat_test(ii, jj));
    }
  }

  Tensor::Matrix9Mandel II = Tensor::constructMatrix9Mandel(I, I);  
  Tensor::Matrix9Mandel II_test; 
  II_test << 1, 1, 1, 0, 0, 0, 0, 0, 0,
             1, 1, 1, 0, 0, 0, 0, 0, 0,
             1, 1, 1, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0;
  for (auto ii = 0u; ii < 9; ++ii) {
    for (auto jj = 0u; jj < 9; ++jj) {
      ASSERT_DOUBLE_EQ(II(ii, jj), II_test(ii, jj));
    }
  }

  Uintah::Vector a(2, 0, -5);
  Uintah::Vector b(8, 3, 2);
  Tensor::Vector9Mandel ab = Tensor::constructVector9Mandel(a, b);  
  double s2 = std::sqrt(2);
  Tensor::Vector9Mandel test_ab;
  test_ab << 16, 0, -10, -15/s2, -18*s2, 3*s2, -15/s2, 22*s2, -3*s2;
  for (auto ii = 0u; ii < 9; ++ii) {
    ASSERT_NEAR(ab(ii), test_ab(ii), 1.0e-6);
  }

  Uintah::Matrix3 C(4, -5, 0, 0, 1, 7, 2, 3, 0);
  Tensor::Vector9Mandel C_vec = Tensor::constructVector9Mandel(C);  
  Tensor::Matrix9Mandel abC = Tensor::constructMatrix9Mandel(a, b, C);  
  Tensor::Matrix9Mandel abC_vec = Tensor::constructMatrix9Mandel(ab, C_vec);  

  //std::cout << "abC = \n" << abC << "\n";
  //std::cout << "abC_vec = \n" << abC_vec << "\n";

  Tensor::Matrix9Mandel test_abC;
  test_abC << 64, 16, 0, 80*s2, 16*s2, -40*s2, -32*s2, -16*s2, 40*s2,
              0, 0, 0, 0, 0, 0, 0, 0, 0,
              -40, -10, 0, -50*s2, -10*s2, 25*s2, 20*s2, 10*s2, -25*s2,
              -30*s2, -15/s2, 0, -75, -15, 75/2.0, 30, 15, -75/2.0,
              -72*s2, -18*s2, 0, -180, -36, 90, 72, 36, -90,
              12*s2, 3*s2, 0, 30, 6, -15, -12, -6, 15,
              -30*s2, -15/s2, 0, -75, -15, 75/2.0, 30, 15, -75/2.0,
              88*s2, 22*s2, 0, 220, 44, -110, -88, -44, 110,
              -12*s2, -3*s2, 0, -30, -6, 15, 12, 6, -15;
  for (auto ii = 0u; ii < 9; ++ii) {
    for (auto jj = 0u; jj < 9; ++jj) {
      ASSERT_NEAR(abC(ii, jj), test_abC(ii, jj), 1.0e-6);
      ASSERT_NEAR(abC_vec(ii, jj), test_abC(ii, jj), 1.0e-6);
    }
  }

  Uintah::Matrix3 ab_mat(a, b);
  Uintah::Matrix3 ab_symm = (ab_mat + ab_mat.Transpose())*0.5;
  Uintah::Matrix3 C_symm = (C + C.Transpose())*0.5;
  Tensor::Matrix9Mandel abC_symm = Tensor::constructMatrix9Mandel(ab_symm, C_symm);  
  Tensor::Matrix6Mandel abC_symm_6 = Tensor::constructMatrix6Mandel(ab_symm, C_symm);  
  //std::cout << "abC_symm = \n" << abC_symm << "\n";
  //std::cout << "abC_symm_6 = \n" << abC_symm_6 << "\n";
  for (auto ii = 0u; ii < 6; ++ii) {
    for (auto jj = 0u; jj < 6; ++jj) {
      ASSERT_NEAR(abC_symm(ii, jj), abC_symm_6(ii, jj), 1.0e-6);
    }
  }

  Tensor::Vector6Mandel C_vec6 = Tensor::constructVector6Mandel(C_symm);
  Uintah::Matrix3 C_back = Tensor::constructMatrix3(C_vec6);
  for (auto ii = 0u; ii < 3; ++ii) {
    for (auto jj = 0u; jj < 3; ++jj) {
      ASSERT_NEAR(C_symm(ii, jj), C_back(ii, jj), 1.0e-6);
    }
  }
  //std::cout << "C_symm = \n" << C_symm << "\n";
  //std::cout << "C_back = \n" << C_back << "\n";

}

TEST(TensorUtilsTest, operations)
{
  Uintah::Matrix3 S(4000, -5000, 1000, -5000, 1000, 7000, 1000, 7000, -5000);
  Uintah::Matrix3 D(1.0e-3, 0, 0,  0, -0.5e-3, 0, 0, 0, -0.5e-3);
  double SD = S.Contract(D);
  Uintah::Matrix3 f(1.0, 2.0, 3.0,  2.0, 4.0, 5.0, 3.0, 5.0, 6.0);
  Tensor::Matrix6Mandel SD_mat = Tensor::constructMatrix6Mandel(S, D);  
  Tensor::Vector6Mandel f_vec = Tensor::constructVector6Mandel(f);  
  Tensor::Matrix6Mandel I = Tensor::IdentityMatrix6Mandel();
  double sig_s = 1000;
  double G = 10000;
  double alpha = sig_s + 2*G*SD;
  auto alphaI = I * alpha;
  auto J = alphaI + 2*G*SD_mat;
  auto J_det = J.determinant();
  auto J_inv = J.inverse();
  auto f_new = J_inv * f_vec;

  ASSERT_DOUBLE_EQ(SD, 6);
  Tensor::Matrix6Mandel SD_check;
  SD_check <<  4, 1, -5, 9.89949, 1.41421, -7.07107,
                  -2, -0.5, 2.5, -4.94975, -0.707107, 3.53553,
                  -2, -0.5, 2.5, -4.94975, -0.707107, 3.53553,
                   0, 0, -0, 0, 0, -0,
                   0, 0, -0, 0, 0, -0,
                   0, 0, -0, 0, 0, -0;
  Tensor::Matrix6Mandel SD_mat_check = SD_check.transpose();

  for (auto ii = 0u; ii < 6; ++ii) {
    for (auto jj = 0u; jj < 6; ++jj) {
      ASSERT_NEAR(SD_mat(ii, jj), SD_mat_check(ii, jj), 1.0e-5);
    }
  }
  ASSERT_NEAR(J_det, 6.25092e+30, 1.0e25);

  Tensor::Vector6Mandel f_check;
  f_check << 1.9238e-05, 3.58012e-05, 3.58698e-05, 8.55968e-05, 3.89429e-05, 3.97673e-06;
  for (auto ii = 0u; ii < 6; ++ii) {
    ASSERT_NEAR(f_new(ii), f_check(ii), 1.0e-5);
  }

  //std::cout << "S = \n" << S << "\n";
  //std::cout << "D = \n" << D << "\n";
  //std::cout << "S:D = " << SD << " alpha = " << alpha << "\n";
  //std::cout << "SD_mat = \n" << SD_mat << "\n";
  //std::cout << "SD_mat_check = \n" << SD_mat_check << "\n";
  //std::cout << "alphaI = \n" << alphaI << "\n";
  //std::cout << "J = \n" << J << "\n";
  //std::cout << "J_det = " << J_det << "\n";
  //std::cout << "J_inv = \n" << J_inv << "\n";
  //std::cout << "f_vec = \n" << f_vec.transpose() << "\n";
  //std::cout << "f_new = \n" << f_new.transpose() << "\n";

  double sigma_m, sigma_s;
  Tensor::Vector6Mandel I_hat, S_hat;
  std::tie(sigma_m, sigma_s, I_hat, S_hat) = Tensor::computeIsomorphicDecomposition(S);

  //std::cout << "I_hat = " << I_hat.transpose() << "\n";
  //std::cout << "s_hat = " << S_hat.transpose() << "\n";
  //std::cout << "sigma_m = " << sigma_m << " sigma_s = " << sigma_s << "\n";

  double I1, sqrt_J2;
  Tensor::Vector6Mandel I_mandel, S_mandel;
  std::tie(I1, sqrt_J2, I_mandel, S_mandel) = Tensor::computeVolDevDecomposition(S);

  //std::cout << "I_mandel = " << I_mandel.transpose() << "\n";
  //std::cout << "s_mandel = " << S_mandel.transpose() << "\n";
  //std::cout << "I1 = " << I1 << "sqrt_J2 = " << sqrt_J2 << "\n";
  
  double sigma_m_check, sigma_s_check;
  double I1_check, sqrt_J2_check;
  Tensor::Vector6Mandel I_hat_check, S_hat_check;
  Tensor::Vector6Mandel I_mandel_check, S_mandel_check;
  I_hat_check << 0.57735, 0.57735, 0.57735, 0, 0, 0;
  S_hat_check << 0.288675, 0.0721688, -0.360844,  0.714435, 0.102062, -0.51031;
  sigma_m_check = 0;
  sigma_s_check = 13856.4;
  I_mandel_check << 1, 1, 1, 0, 0, 0;
  S_mandel_check << 4000, 1000, -5000, 9899.49, 1414.21, -7071.07;
  I1_check = 0;
  sqrt_J2_check = 9797.96;

  for (auto ii = 0u; ii < 6; ++ii) {
    ASSERT_NEAR(I_hat(ii), I_hat_check(ii), 1.0e-5);
    ASSERT_NEAR(S_hat(ii), S_hat_check(ii), 1.0e-5);
    ASSERT_NEAR(I_mandel(ii), I_mandel_check(ii), 1.0e-5);
    ASSERT_NEAR(S_mandel(ii), S_mandel_check(ii), 1.0e-1);
    ASSERT_NEAR(S_hat(ii), S_hat_check(ii), 1.0e-5);
  }
  ASSERT_NEAR(sigma_m, sigma_m_check, 1.0e-5);
  ASSERT_NEAR(sigma_s, sigma_s_check, 1.0e-2);
  ASSERT_NEAR(I1, I1_check, 1.0e-5);
  ASSERT_NEAR(sqrt_J2, sqrt_J2_check, 1.0e-2);
  
}
