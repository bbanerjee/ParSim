#include <CCA/Components/MPM/ConstitutiveModel/Models/TensorUtils.h>

#include <iostream>
#include <vector>

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

}

