#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Default.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arena.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <typeinfo>

#include <gtest/gtest.h>

struct TestPressureModel 
{
  double
  computePressureDynamic(const Vaango::ModelStateBase* state_in)
  {
    auto state = dynamic_cast<const Vaango::ModelState_Default*>(state_in);
    if (!state) {
      std::cout << "Dynamic cast: Wrong state object\n";
      return 0.0;
    }
    return state->pressure;
  }

  double
  computePressureStatic(const Vaango::ModelStateBase* state_in)
  {
    auto state = static_cast<const Vaango::ModelState_Default*>(state_in);
    return state->pressure;
  }

  double
  computePressureTypeid(const Vaango::ModelStateBase* state_in)
  {
    double pressure = 0.0;
    if (typeid(state_in) != typeid(Vaango::ModelState_Default*)) {
      //std::cout << "Typeid: Wrong state object\n";
      //std::cout << "state_in has type: " << typeid(state_in).name() << "\n";
      //std::cout << "compare has type: " << typeid(Vaango::ModelState_Default*).name() << "\n";
      pressure = 1.0; 
    }
    pressure += 1.0;

    auto state = static_cast<const Vaango::ModelState_Default*>(state_in);
    return state->pressure + pressure;
  }

  double
  computePressureNumVar(const Vaango::ModelStateBase* state_in)
  {
    double pressure = 0.0;
    Vaango::ModelState_Default s;
    if (state_in->numStateVar() != s.numStateVar()) {
      //std::cout << "NumStateVar: Wrong state object\n";
      //std::cout << "state_in has vars: " << state_in->numStateVar() << "\n";
      //std::cout << "compare has vars: " << s.numStateVar() << "\n";
      pressure = 1.0; 
    }

    auto state = static_cast<const Vaango::ModelState_Default*>(state_in);
    return state->pressure;
  }
};

TEST(ModelStateTest, Casts)
{
  Vaango::ModelState_Default state;
  state.bulkModulus = 1.0e6;
  state.pressure = 1.0e3;

  TestPressureModel p_model;

  double p_dynamic = 0.0;
  double p_static = 0.0;
  double p_typeid = 0.0;
  double p_numvar = 0.0;

  auto start_dynamic = std::chrono::steady_clock::now();
  int count = 0;
  while (count < 10000000) {
    p_dynamic = p_model.computePressureDynamic(&state);
    ++count;
  }
  auto end_dynamic = std::chrono::steady_clock::now();
  
  auto start_static = std::chrono::steady_clock::now();
  count = 0;
  while (count < 10000000) {
    p_static = p_model.computePressureStatic(&state);
    ++count;
  }
  auto end_static = std::chrono::steady_clock::now();

  auto start_typeid = std::chrono::steady_clock::now();
  count = 0;
  while (count < 10000000) {
    p_typeid = p_model.computePressureTypeid(&state);
    ++count;
  }
  auto end_typeid = std::chrono::steady_clock::now();

  auto start_numvar = std::chrono::steady_clock::now();
  count = 0;
  while (count < 10000000) {
    p_numvar = p_model.computePressureNumVar(&state);
    ++count;
  }
  auto end_numvar = std::chrono::steady_clock::now();

  std::chrono::duration<double> dynamic_seconds = end_dynamic - start_dynamic;
  std::chrono::duration<double> static_seconds = end_static - start_static;
  std::chrono::duration<double> typeid_seconds = end_typeid - start_typeid;
  std::chrono::duration<double> numvar_seconds = end_numvar - start_numvar;
  std::cout << "Dynamic: elapsed time: " << dynamic_seconds.count() << "s\n";
  std::cout << "Static: elapsed time: " << static_seconds.count() << "s\n";
  std::cout << "TypeID: elapsed time: " << typeid_seconds.count() << "s\n";
  std::cout << "NumVar: elapsed time: " << numvar_seconds.count() << "s\n";

  ASSERT_DOUBLE_EQ(p_dynamic, state.pressure);
  ASSERT_DOUBLE_EQ(p_static, state.pressure);
  ASSERT_DOUBLE_EQ(p_typeid, state.pressure + 2);
  // ASSERT_NEAR(errMat(ii, jj), 0.0, 1.0e-10);

  Vaango::ModelState_Arena state_arena;
  double p_dynamic_arena = p_model.computePressureDynamic(&state_arena);
  ASSERT_DOUBLE_EQ(p_dynamic_arena, state.pressure * 0.0);
  double p_static_arena = p_model.computePressureStatic(&state_arena);
  ASSERT_DOUBLE_EQ(p_static_arena, state.pressure * 0.0);
  double p_numvar_arena = p_model.computePressureNumVar(&state_arena);
  ASSERT_DOUBLE_EQ(p_numvar_arena, state.pressure * 0.0);
}

