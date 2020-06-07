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
  ASSERT_NEAR(denominator, 22466, 1.0);

  Vector6 dEps;
  dEps  << 1, 2, 3, 4, 5, 6;
  //std::cout << "dEps = \n" << dEps << "\n";

  double numerator_new = df_dsigma_T * dEps;
  //std::cout << "numerator_new = \n" << numerator_new << "\n";
  ASSERT_NEAR(numerator_new, 53660, 1.0);

  double numerator_other = numerator * dEps;
  //std::cout << "numerator_other = \n" << numerator_other << "\n";
  ASSERT_NEAR(numerator_new, numerator_other, 1.0);

  auto num = df_dot_D * dEps;
  auto den = df_dot_D * df_dsigma;
  //std::cout << "num = " << num << " den = " << den << "\n";
  ASSERT_NEAR(num, numerator_new, 1.0);
  ASSERT_NEAR(den, denominator, 1.0);

  
  Vector6 stress;
  stress << 1100, 2200, 3300, 700, 500, 600;

  using Vector3 = Eigen::Matrix<double, 3, 1>;
  using Matrix33 = Eigen::Matrix<double, 3, 3>;
  Matrix33 stressMat;
  stressMat(0, 0) = stress(0);
  stressMat(0, 1) = stress(3);
  stressMat(0, 2) = stress(4);
  stressMat(1, 0) = stressMat(0, 1); 
  stressMat(1, 1) = stress(1);
  stressMat(1, 2) = stress(5);
  stressMat(2, 0) = stressMat(0, 2); 
  stressMat(2, 1) = stressMat(1, 2);
  stressMat(2, 2) = stress(2);
  
  Eigen::SelfAdjointEigenSolver<Matrix33> solver(stressMat);
  Vector3 eigenval = solver.eigenvalues();
  Matrix33 eigenvec = solver.eigenvectors();

  //std::cout << "eigenvals = " << eigenval << "\n";
  //std::cout << "eigenvecs = " << eigenvec << "\n";
  
  double temp_val = eigenval(0);
  eigenval(0) = eigenval(2);
  eigenval(2) = temp_val;
  //std::cout << "eigenvals = " << eigenval << "\n";

  auto temp_vec = eigenvec.col(0).eval();
  eigenvec.col(0) = eigenvec.col(2);
  eigenvec.col(2) = temp_vec;
  //std::cout << "eigenvecs = " << eigenvec << "\n";

  Vector6 principal = Vector6::Zero();
  principal(0) = eigenval(0);
  principal(1) = eigenval(1);
  principal(2) = eigenval(2);
  //std::cout << "principal = " << principal << "\n";

  Matrix33 prinMat;
  prinMat(0, 0) = principal(0);
  prinMat(0, 1) = principal(3);
  prinMat(0, 2) = principal(4);
  prinMat(1, 0) = prinMat(0, 1); 
  prinMat(1, 1) = principal(1);
  prinMat(1, 2) = principal(5);
  prinMat(2, 0) = prinMat(0, 2); 
  prinMat(2, 1) = prinMat(1, 2);
  prinMat(2, 2) = principal(2);
  //std::cout << "prinMat = " << prinMat << "\n";

  auto rotMat = eigenvec * prinMat * eigenvec.transpose();
  //std::cout << rotMat << "\n";

  auto errMat = rotMat - stressMat;
    
  for (auto ii = 0u; ii < 3; ++ii) {
    for (auto jj = 0u; jj < 3; ++jj) {
      ASSERT_NEAR(errMat(ii, jj), 0.0, 1.0e-10);
    }
  }

  //Vector3 one = Vector3::Ones();
  //std::cout << one << "\n";

  Vector3 a, b;
  a << 1, 2, 3;
  b << 4, 5, 6;

  auto axb = a.cross(b);
  //std::cout << "axb = " << axb.transpose() << "\n";

  Vector3 cross_prod;
  cross_prod << a(1) * b(2) - b(1) * a(2),
                -a(0) * b(2) + b(0) * a(2),
                a(0) * b(1) - b(0) * a(1);
  //std::cout << "cross_prod = " << cross_prod.transpose() << "\n";

  for (auto ii = 0u; ii < 3; ++ii) {
    ASSERT_DOUBLE_EQ(axb(ii), cross_prod(ii));
  }
  
  /* Test abs and max */
  Vector7 strain_final, strain_initial;
  strain_final << 0.0129184, -8.83862e-39, -2.96308e-38, 
                 -1.36816e-19, -2.37293e-19, -1.4483e-38, 0;
  strain_initial << 0.012923, -8.24534e-39, -2.99804e-38, 
                    -1.30285e-19, -2.39469e-19, -1.40575e-38, 0;
  Vector7 strain_inc = strain_final - strain_initial;
  double max_strain = strain_inc.block<6,1>(0,0).cwiseAbs().maxCoeff();
  //std::cout << max_strain <<"\n";

  double max = 0.0;
  for (int i = 0; i < 6; i++) {
    if (std::abs(strain_inc(i)) > max) {
      max = std::abs(strain_inc(i));
    }
  }
  //std::cout << max << "\n";
  ASSERT_DOUBLE_EQ(max_strain, max);

  Vector6 eps_inc_1 = strain_inc.block<6,1>(0,0) / (max_strain * 1.0e-10);
  Vector6 eps_inc_2 = strain_inc.block<6,1>(0,0) / (max * 1.0e-10);
  
  ASSERT_DOUBLE_EQ(eps_inc_1.squaredNorm(), eps_inc_2.squaredNorm());
  //std::cout << eps_inc_1.transpose() << "\n"
  //          << eps_inc_2.transpose() << "\n";
}

