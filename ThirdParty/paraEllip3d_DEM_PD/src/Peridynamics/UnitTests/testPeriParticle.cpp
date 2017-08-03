#include <Peridynamics/PeriParticle.h>
#include <Core/Util/Utility.h>
#include <InputOutput/InputParameter.h>
#include <gtest/gtest.h>

using namespace pd;

TEST(PeriParticleTest, accFunctions) {

  // Setup the paramters that are used by the constructor
  dem::InputParameter::get().addParameter("typeConstitutive", 1.0);
  dem::InputParameter::get().addParameter("chi", 1.0);
  dem::InputParameter::get().addParameter("tangentModulus11", 1.0);
  dem::InputParameter::get().addParameter("tangentModulus12", 1.0);
  dem::InputParameter::get().addParameter("tangentModulus13", 1.0);
  dem::InputParameter::get().addParameter("tangentModulus21", 1.0);
  dem::InputParameter::get().addParameter("tangentModulus22", 1.0);
  dem::InputParameter::get().addParameter("tangentModulus23", 1.0);
  dem::InputParameter::get().addParameter("tangentModulus31", 1.0);
  dem::InputParameter::get().addParameter("tangentModulus32", 1.0);
  dem::InputParameter::get().addParameter("tangentModulus33", 1.0);
  dem::InputParameter::get().addParameter("tangentModulus44", 1.0);
  dem::InputParameter::get().addParameter("tangentModulus55", 1.0);
  dem::InputParameter::get().addParameter("tangentModulus66", 1.0);
  PeriParticle particle;
  EXPECT_EQ(particle.getId(), 0);
}
