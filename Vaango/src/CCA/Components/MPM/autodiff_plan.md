# Implementation Plan: Automatic Differentiation for Consistent Tangent Stiffness

This document outlines the strategy for incorporating **CoDiPack** to compute consistent material tangent stiffnesses ($D_{ijkl} = \frac{\partial \sigma_{ij}}{\partial \epsilon_{kl}}$ or $\frac{\partial P_{iJ}}{\partial F_{kL}}$) in Vaango's implicit constitutive models.

## 1. Build System Integration
- **Dependency**: Add [CoDiPack](https://github.com/SciCompInst/CoDiPack) as a git submodule in `src/submodules/codipack`.
- **CMake**: Update the root `CMakeLists.txt` and `CCA/Components/MPM/CMakeLists.txt` to include the CoDiPack header directory.
    - Note: CoDiPack is header-only, so no library linking is required.
    - Ensure C++17 or higher is enabled (required by modern AD libraries).

## 2. Refactoring Constitutive Models
To use AD, the stress computation must be generic with respect to the scalar type.

### Step A: Templating Stress Logic
Refactor the core stress calculation into a templated function (or functor) in the constitutive model class.
```cpp
template <typename T>
void computeStressTemplated(const Eigen::Matrix<T, 3, 3>& F, 
                            const MaterialParams<T>& params,
                            Eigen::Matrix<T, 3, 3>& stress);
```

### Step B: The AD Wrapper
In the standard `computeTangentStiffness` (which uses `double`), implement the CoDiPack forward-mode routine:
1. **Initialize Active Types**: Use `codi::RealForward` as the scalar type for inputs.
2. **Seeding**: Seed the deformation gradient (or strain) components.
3. **Execution**: Call `computeStressTemplated`.
4. **Extraction**: Retrieve the partial derivatives (Jacobian) from the active output stress.

## 3. Mapping to Voigt Form (Matrix66)
MPM implicit solvers typically expect a 6x6 tangent modulus ($D$). 
- Map the derivatives of the 3x3 stress tensor with respect to the 3x3 deformation/strain tensor into the symmetric 6x6 Voigt representation.
- Ensure the mapping accounts for the symmetry of the Cauchy stress ($\sigma_{ij} = \sigma_{ji}$) and the specific Voigt ordering used in Vaango (usually 11, 22, 33, 23, 13, 12).

## 4. Implicit Solver Interaction (`ImpMPM.cc`)
- The current `formStiffnessMatrix` and `computeInternalForce` in `ImpMPM.cc` call `computeStressTensorImplicit`.
- The `ImplicitCM` interface already passes a `Solver*` pointer. The AD-derived tangents will be filled into the global stiffness matrix using the existing `BtDB` and `BnltDBnl` mechanisms in `ImplicitCM.cc`.

## 5. Validation Strategy
1. **Unit Tests**: Create a test case in `CCA/Components/MPM/ConstitutiveModel/Utilities/UnitTests` to compare AD-derived tangents with analytical tangents for a simple model (e.g., `CompNeoHook`).
2. **Convergence Rates**: Monitor the residual convergence in `ImpMPM`. A truly consistent tangent should yield quadratic convergence in the Newton-Raphson iterations.
3. Regression**: Ensure `double`-based stress results remain identical to the previous implementation.

---

# Alternative Plan: Using the 'autodiff' Library

This is a modern C++17 alternative that provides deeper integration with **Eigen**, which is already heavily used in the implicit MPM components.

## 1. Library Overview
- **Library**: [autodiff](https://autodiff.github.io/) (Forward-mode)
- **Strengths**: 
    - Native support for Eigen types (`Eigen::Matrix<autodiff::real, ...>`).
    - Very clean high-level API for Jacobians and Hessians.
    - Header-only, minimal configuration.

## 2. Implementation Details

### A. Integrated Headers
Include the library in your constitutive model:
```cpp
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
```

### B. Templated Physics Core
Maintain a single templated function for the constitutive logic:
```cpp
template <typename T>
Eigen::Matrix<T, 3, 3> 
computeCauchyStress(const Eigen::Matrix<T, 3, 3>& F, ...) {
    // This logic works for both 'double' and 'autodiff::real'
}
```

### C. Computing the Material Tangent
Use the `autodiff::jacobian` function to differentiate the stress with respect to the deformation gradient $F$ or strain $E$.

```cpp
Matrix66 computeTangentStiffness(const Matrix3& F_uintah) {
    // 1. Convert to Eigen
    Eigen::Matrix3d F_val = convertToEigen(F_uintah);

    // 2. Define functional wrapper for differentiation
    auto stress_func = [&](const Eigen::VectorXd& F_vec) -> Eigen::VectorXd {
        auto F_active = F_vec.reshaped(3, 3).cast<autodiff::real>();
        auto sigma = computeCauchyStress(F_active, ...);
        return mapToVoigt(sigma);
    };

    // 3. Automated Jacobian calculation
    Eigen::VectorXd F_flat = F_val.reshaped();
    Eigen::VectorXd sigma_voigt;
    Eigen::MatrixXd J = autodiff::jacobian(stress_func, wrt(F_flat), at(F_flat), sigma_voigt);

    // 4. J is the consistent tangent (9x9 or 6x9), map to Matrix66
    return reduceTo6x6(J);
}
```

## 3. Comparison with CoDiPack
- **Syntax**: `autodiff` is more "functional" and feels like standard C++, while CoDiPack is slightly more imperative (seeding/unseeding).
- **Performance**: CoDiPack is generally faster for very large-scale reverse-mode problems, but for the small-scale 9x9 Jacobians required for material tangents, the difference is negligible.
- **Ease of Use**: `autodiff`'s native Eigen support makes it the "easiest" path for immediate implementation in Vaango's existing Eigen-based models.

