#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Tabular.h>

#include <iostream>

#include <gtest/gtest.h>

using namespace Vaango;

TEST(ModelStateTabularTest, DefaultConstructor)
{
  ModelState_Tabular state;
  EXPECT_EQ(state.particleID, 0);
  EXPECT_DOUBLE_EQ(state.stressTensor(0,0), 0);
  EXPECT_DOUBLE_EQ(state.deviatoricStressTensor(0,0), 0);
  EXPECT_DOUBLE_EQ(state.elasticStrainTensor(0,0), 0);
  EXPECT_DOUBLE_EQ(state.plasticStrainTensor(0,0), 0);
  EXPECT_DOUBLE_EQ(state.bulkModulus, 0);
  EXPECT_DOUBLE_EQ(state.shearModulus, 0);
}

TEST(ModelStateTabularTest, CopyConstructor)
{
  ModelState_Tabular state;
  state.particleID = 1234567890123;
  state.stressTensor = Uintah::Matrix3(1000, 2000, 3000, 
                                       2000, 4000, 5000,
                                       3000, 5000, 6000);
  state.updateStressInvariants();
  state.elasticStrainTensor = Uintah::Matrix3(0.1, 0.2, 0.3,
                                              0.2, 1.0, 2.0,
                                              0.3, 2.0, 5.0);
  state.plasticStrainTensor = Uintah::Matrix3(0.2, 0.3, 0.4,
                                              0.3, 2.0, 3.0,
                                              0.4, 3.0, 6.0);
  state.updatePlasticStrainInvariants();
  EXPECT_NEAR(state.deviatoricStressTensor.Trace(), 0, 1.0e-12);

  ModelState_Tabular stateCopy(state);
  EXPECT_DOUBLE_EQ(state.stressTensor(0,0), stateCopy.stressTensor(0,0));
  EXPECT_DOUBLE_EQ(state.deviatoricStressTensor(0,0), stateCopy.deviatoricStressTensor(0,0));
  EXPECT_DOUBLE_EQ(state.elasticStrainTensor(0,0), stateCopy.elasticStrainTensor(0,0));
  EXPECT_DOUBLE_EQ(state.plasticStrainTensor(0,0), stateCopy.plasticStrainTensor(0,0));
  EXPECT_DOUBLE_EQ(state.bulkModulus, stateCopy.bulkModulus);
  EXPECT_DOUBLE_EQ(state.shearModulus, stateCopy.shearModulus);
}

TEST(ModelStateTabularTest, MoveConstructor)
{
  ModelState_Tabular state;
  state.particleID = 1234567890123;
  state.stressTensor = Uintah::Matrix3(1000, 2000, 3000, 
                                       2000, 4000, 5000,
                                       3000, 5000, 6000);
  state.updateStressInvariants();
  state.elasticStrainTensor = Uintah::Matrix3(0.1, 0.2, 0.3,
                                              0.2, 1.0, 2.0,
                                              0.3, 2.0, 5.0);
  state.plasticStrainTensor = Uintah::Matrix3(0.2, 0.3, 0.4,
                                              0.3, 2.0, 3.0,
                                              0.4, 3.0, 6.0);
  state.updatePlasticStrainInvariants();
  EXPECT_NEAR(state.deviatoricStressTensor.Trace(), 0, 1.0e-12);

  ModelState_Tabular stateMove(std::move(state));
  EXPECT_DOUBLE_EQ(state.stressTensor(0,0), stateMove.stressTensor(0,0));
  EXPECT_DOUBLE_EQ(state.deviatoricStressTensor(0,0), stateMove.deviatoricStressTensor(0,0));
  EXPECT_DOUBLE_EQ(state.elasticStrainTensor(0,0), stateMove.elasticStrainTensor(0,0));
  EXPECT_DOUBLE_EQ(state.plasticStrainTensor(0,0), stateMove.plasticStrainTensor(0,0));
  EXPECT_DOUBLE_EQ(state.bulkModulus, stateMove.bulkModulus);
  EXPECT_DOUBLE_EQ(state.shearModulus, stateMove.shearModulus);
}

TEST(ModelStateTabularTest, Assignment)
{
  ModelState_Tabular state;
  state.particleID = 1234567890123;
  state.stressTensor = Uintah::Matrix3(1000, 2000, 3000, 
                                       2000, 4000, 5000,
                                       3000, 5000, 6000);
  state.updateStressInvariants();
  state.elasticStrainTensor = Uintah::Matrix3(0.1, 0.2, 0.3,
                                              0.2, 1.0, 2.0,
                                              0.3, 2.0, 5.0);
  state.plasticStrainTensor = Uintah::Matrix3(0.2, 0.3, 0.4,
                                              0.3, 2.0, 3.0,
                                              0.4, 3.0, 6.0);
  state.updatePlasticStrainInvariants();
  EXPECT_NEAR(state.deviatoricStressTensor.Trace(), 0, 1.0e-12);

  ModelState_Tabular stateCopy = state;
  EXPECT_DOUBLE_EQ(state.stressTensor(0,0), stateCopy.stressTensor(0,0));
  EXPECT_DOUBLE_EQ(state.deviatoricStressTensor(0,0), stateCopy.deviatoricStressTensor(0,0));
  EXPECT_NEAR(stateCopy.deviatoricStressTensor.Trace(), 0, 1.0e-12);
  EXPECT_DOUBLE_EQ(state.elasticStrainTensor(0,0), stateCopy.elasticStrainTensor(0,0));
  EXPECT_DOUBLE_EQ(state.plasticStrainTensor(0,0), stateCopy.plasticStrainTensor(0,0));
  EXPECT_DOUBLE_EQ(state.bulkModulus, stateCopy.bulkModulus);
  EXPECT_DOUBLE_EQ(state.shearModulus, stateCopy.shearModulus);
}

TEST(ModelStateTabularTest, AssignmentPointer)
{
  ModelState_Tabular state;
  state.particleID = 1234567890123;
  state.stressTensor = Uintah::Matrix3(1000, 2000, 3000, 
                                       2000, 4000, 5000,
                                       3000, 5000, 6000);
  state.updateStressInvariants();
  state.elasticStrainTensor = Uintah::Matrix3(0.1, 0.2, 0.3,
                                              0.2, 1.0, 2.0,
                                              0.3, 2.0, 5.0);
  state.plasticStrainTensor = Uintah::Matrix3(0.2, 0.3, 0.4,
                                              0.3, 2.0, 3.0,
                                              0.4, 3.0, 6.0);
  state.updatePlasticStrainInvariants();
  EXPECT_NEAR(state.deviatoricStressTensor.Trace(), 0, 1.0e-12);

  ModelState_Tabular* stateCopy = &state;
  //std::cout << *stateCopy;
  EXPECT_DOUBLE_EQ(state.stressTensor(0,0), stateCopy->stressTensor(0,0));
  EXPECT_DOUBLE_EQ(state.deviatoricStressTensor(0,0), stateCopy->deviatoricStressTensor(0,0));
  EXPECT_NEAR(stateCopy->deviatoricStressTensor.Trace(), 0, 1.0e-12);
  EXPECT_DOUBLE_EQ(state.elasticStrainTensor(0,0), stateCopy->elasticStrainTensor(0,0));
  EXPECT_DOUBLE_EQ(state.plasticStrainTensor(0,0), stateCopy->plasticStrainTensor(0,0));
  EXPECT_DOUBLE_EQ(state.bulkModulus, stateCopy->bulkModulus);
  EXPECT_DOUBLE_EQ(state.shearModulus, stateCopy->shearModulus);
}
