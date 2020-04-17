#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include <gtest/gtest.h>

TEST(EigenTest, constructors)
{
  using Vector7 = Eigen::Matrix<double, 7, 1, Eigen::DontAlign>;
  Vector7 reuseRes = Vector7::Zero();

  Eigen::Matrix<double, 6, 2> plasticStrain = Eigen::Matrix<double, 6, 2>::Zero();

  //std::cout << "after init:" << reuseRes << "\n";
  //std::cout << "after init:" << plasticStrain << "\n";

  reuseRes << 101, 102, 103, 104, 105, 106, 107;
  plasticStrain << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;

  //std::cout << "after set:" << reuseRes << "\n";
  //std::cout << "after set:" << plasticStrain << "\n";

  plasticStrain.col(0) = reuseRes.block<6,1>(0,0);
  //std::cout << "after block op:" << plasticStrain << "\n";

  //std::cout << "last elem = " << reuseRes(6) << "\n";

  double K = 1000;
  double G = 700;
  double K43G = K + 4.0 * G / 3.0;
  double K23G = K - 2.0 * G / 3.0;

  using Matrix66 = Eigen::Matrix<double, 6, 6, Eigen::DontAlign>;
  Matrix66 DEP = Matrix66::Zero();
  DEP(0, 0) =  K43G;
  DEP(0, 1) =  K23G;
  DEP(0, 2) =  K23G;
  DEP(1, 0) =  K23G;
  DEP(1, 1) =  K43G;
  DEP(1, 2) =  K23G;
  DEP(2, 0) =  K23G;
  DEP(2, 1) =  K23G;
  DEP(2, 2) =  K43G;
  DEP(3, 3) =  2.0 * G;
  DEP(4, 4) =  2.0 * G;
  DEP(5, 5) =  2.0 * G;
  //std::cout << "DEP = \n" << DEP << "\n";

  using Matrix67 = Eigen::Matrix<double, 6, 7>;
  Matrix67 DEP_extra = Matrix67::Zero();
  DEP_extra.block<6,6>(0,0) = DEP;
  //std::cout << "DEP_extra = \n" << DEP_extra << "\n";
  
  using Vector6 = Eigen::Matrix<double, 6, 1>;
  Vector6 df_dsigma;
  df_dsigma << 1.100, 1.200, 1.300, 1.400, 1.500, 1.600;

  using Vector6T = Eigen::Matrix<double, 1, 6>;
  Vector6T df_dsigma_T = df_dsigma.transpose();
  //std::cout << "df_dsigma_T = \n" << df_dsigma_T << "\n";
  
  Vector6T numerator = df_dsigma_T * DEP;
  //std::cout << "numerator = \n" << numerator << "\n";

  Vector6T df_dot_D = df_dsigma.transpose() * DEP;
  //std::cout << "df_dot_D = \n" << df_dot_D << "\n";

  df_dsigma_T *= DEP;
  //std::cout << "df_dsigma_T = \n" << df_dsigma_T << "\n";

  double denominator = numerator * df_dsigma;
  //std::cout << "denominator = \n" << denominator << "\n";

  Vector6 dEps;
  dEps  << 1, 2, 3, 4, 5, 6;
  //std::cout << "dEps = \n" << dEps << "\n";

  double numerator_new = df_dsigma_T * dEps;
  //std::cout << "numerator_new = \n" << numerator_new << "\n";

  double numerator_other = numerator * dEps;
  //std::cout << "numerator_other = \n" << numerator_other << "\n";

  auto num = df_dot_D * dEps;
  auto den = df_dot_D * df_dsigma;
  //std::cout << "num = " << num << " den = " << den << "\n";

  
  /*
  for (auto ii = 0u; ii < 6; ++ii) {

    //ASSERT_DOUBLE_EQ(I_vec(ii), I_vec_test(ii));
  }

  for (auto ii = 0u; ii < 9; ++ii) {
    //ASSERT_NEAR(ab(ii), test_ab(ii), 1.0e-6);
  }
  */

}

