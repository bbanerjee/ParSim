# Implementation Plan: Autodiff for Mohr-Coulomb Model

This document describes the steps required to enable implicit support and automatic differentiation for the Mohr-Coulomb constitutive model in Vaango.

## 1. Prerequisites
- **Implicit Support**: The `MohrCoulomb` class currently throws an exception for implicit tasks. `addComputesAndRequires` and `computeStressTensorImplicit` must be implemented.
- **Library**: The `autodiff` header-only library must be available in the include path.

## 2. Templating the Core Logic
To allow `autodiff` to track derivatives through the complex integration algorithms (Runge-Kutta, etc.), the following must be templated:

### A. Data Types (`MohrCoulombTypes.h`)
- Convert `Vector6`, `Vector7`, `Matrix66`, etc., to templated aliases or use `Eigen::Matrix<T, ...>`.
```cpp
template <typename T>
using Vector6T = Eigen::Matrix<T, 6, 1>;
```

### B. State Object (`MohrCoulombState.h`)
- Template `MohrCoulombState` so it can hold active types during the tangent calculation.
```cpp
template <typename T>
class MohrCoulombStateT {
  Vector6T<T> stress;
  Vector7T<T> strain;
  // ...
};
```

### C. Base Class and Integration (`MohrCoulombBase.h`)
- Template the `integrate` method and all protected helper methods (`plasticRK543`, `computeYieldNormalized`, etc.).
- This is the most significant change, as it allows the entire integration loop to be differentiated.

## 3. Consistent Tangent Calculation
In `MohrCoulomb::computeStressTensorImplicit`:

1. **Setup**: Create a lambda or functor that wraps the `integrate` call.
2. **Jacobian**: Use `autodiff::jacobian` to compute:
   $$\mathbf{D}_{consistent} = \frac{\partial \Delta \boldsymbol{\sigma}}{\partial \Delta \boldsymbol{\epsilon}}$$
3. **Seeding**: The strain increment $\Delta \boldsymbol{\epsilon}$ is the independent variable.
4. **Result**: The output is a 6x6 (or 6x7 if suction is included) matrix representing the consistent tangent modulus.

## 4. Integration with Implicit Solver
- The derived `D` matrix is passed to the `Solver*` via `BtDB` and `BnltDBnl` in `ImplicitCM`.
- Because `autodiff` differentiates through the *actual* integration path taken (including any sub-stepping or drift correction), the resulting tangent is perfectly consistent with the residual, ensuring quadratic convergence for Newton iterations.

## 5. Implementation Sequence
1.  **Refactor MohrCoulombTypes.h** to support templated Eigen types.
2.  **Refactor MohrCoulombState.h** to `MohrCoulombStateT<T>`.
3.  **Template MohrCoulombBase** and ensure all mathematical operations (sin, cos, sqrt) use `std::` or `autodiff::` overloads.
4.  **Implement implicit hooks** in `MohrCoulomb.cc`.
5.  **Replace manual `getTangentMatrix` calls** with the AD version in the implicit loop.

---

# Alternative Plan: Using the 'CoDiPack' Library

This plan provides an alternative using **CoDiPack** for the Mohr-Coulomb model, which is highly efficient for complex integration algorithms with many intermediate states.

## 1. Prerequisites
- **Library**: CoDiPack should be included as a submodule in `src/submodules/codipack`.
- **Active Type**: Use `codi::RealForward` for forward-mode differentiation.

## 2. Structural Changes
The templating requirements are identical to the `autodiff` plan (see Section 2 above). `MohrCoulombStateT<T>` and `MohrCoulombBase` must be templated to allow `T = codi::RealForward`.

## 3. Forward-Mode AD Routine
In `MohrCoulomb::computeStressTensorImplicit`, the consistent tangent is calculated as follows:

1.  **Input Seeding**:
    *   Initialize the input strain increment $\Delta \boldsymbol{\epsilon}$ as a `Vector7T<codi::RealForward>`.
    *   For each component $i \in [0, 6]$, set the value and seed the derivative:
        ```cpp
        delta_epsilon[i].setValue(double_val);
        delta_epsilon[i].setGradient(1.0); // Seed for component i
        ```
2.  **Functional Evaluation**:
    *   Call the templated `integrate` method with the active types.
    *   Since CoDiPack forward-mode typically differentiates one input at a time (or use a vector mode), you may loop over the 7 components of the strain increment to build the 6x7 Jacobian.

3.  **Gradient Harvesting**:
    *   After each integration, extract the derivatives from the output stress $\boldsymbol{\sigma}$:
        ```cpp
        tangent_column[i] = stress[j].getGradient();
        ```
    *   Reset gradients before the next input component is seeded.

## 4. Why CoDiPack for Mohr-Coulomb?
- **Robustness**: Mohr-Coulomb integration involves complex branching logic (yield surface checks, sub-stepping, different algorithms). CoDiPack is specifically designed to handle such "active" branching correctly.
- **Performance**: For the specific sub-stepping algorithms used in `MohrCoulombBase`, CoDiPack's low-overhead active types provide excellent performance during the tangent matrix assembly.
- **Flexibility**: If the model is later extended to include more state variables, CoDiPack can easily scale to compute those derivatives.

