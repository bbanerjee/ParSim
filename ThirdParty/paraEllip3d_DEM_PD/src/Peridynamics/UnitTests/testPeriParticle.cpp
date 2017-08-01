#include <Peridynamics/PeriParticle.h>
#include <Core/Util/Utility.h>
#include <InputOutput/Parameter.h>
#include <gtest/gtest.h>

using namespace pd;

TEST(PeriParticleTest, accFunctions) {

  // Setup the paramters that are used by the constructor
  dem::Parameter::get().addParameter("typeConstitutive", 1.0);
  dem::Parameter::get().addParameter("chi", 1.0);
  dem::Parameter::get().addParameter("tangentModulus11", 1.0);
  dem::Parameter::get().addParameter("tangentModulus12", 1.0);
  dem::Parameter::get().addParameter("tangentModulus13", 1.0);
  dem::Parameter::get().addParameter("tangentModulus21", 1.0);
  dem::Parameter::get().addParameter("tangentModulus22", 1.0);
  dem::Parameter::get().addParameter("tangentModulus23", 1.0);
  dem::Parameter::get().addParameter("tangentModulus31", 1.0);
  dem::Parameter::get().addParameter("tangentModulus32", 1.0);
  dem::Parameter::get().addParameter("tangentModulus33", 1.0);
  dem::Parameter::get().addParameter("tangentModulus44", 1.0);
  dem::Parameter::get().addParameter("tangentModulus55", 1.0);
  dem::Parameter::get().addParameter("tangentModulus66", 1.0);
  PeriParticle particle;
  EXPECT_EQ(particle.getId(), 0);
}
