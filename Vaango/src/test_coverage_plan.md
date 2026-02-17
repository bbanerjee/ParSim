# Test Coverage Enhancement Plan for Vaango

This document outlines the proposed strategy for increasing the test coverage of the Vaango project, ensuring better code quality and reliability for its multi-physics simulation capabilities.

## Phase 1: Baseline Measurement
Before adding new tests, we must establish a clear understanding of current coverage levels across all modules.
*   **Unified Reporting**: Run both unit tests (compile-time) and regression tests (`R_Tester`) with coverage instrumentation enabled.
*   **Identify Gaps**: Generate a combined HTML report to pinpoint "dark" areas—components with low or zero coverage (e.g., specific ICE solvers or complex MPMICE coupling logic).

## Phase 2: Strengthening the Foundation (Core)
The reliability of simulations depends on the correctness of low-level utilities.
*   **Core/Math**: Ensure 100% coverage for tensor operations, rotation matrices, and interpolation functions. These are used in almost every computational step.
*   **Core/Grid**: Expand tests for patch navigation, boundary condition application, and grid variable management.
*   **Refactoring for Testability**: Identify large, monolithic functions in `Core` and break them into smaller, "pure" functions that can be unit-tested without complex state setup.

## Phase 3: Component-Level Testing (CCA)
Focus on the physics and solver logic.
*   **Constitutive Model Suite**:
    *   Implement a **Parameterised Test Framework** using GTest to run standard stress-strain tests (e.g., uniaxial tension, pure shear) across all material models.
    *   Verify thermodynamic requirements (e.g., non-negative plastic dissipation).
*   **Mocking Infrastructure**:
    *   Introduce `gmock` to simulate components like `DataWarehouse` or `SimulationState`. This allows testing solver logic in isolation from the full parallel runtime.
*   **Boundary Conditions**: Create targeted tests for complex BCs (e.g., symmetry, periodic, custom pressure profiles) which are often only tested via large integration runs.

## Phase 4: Integration and Validation (Regression)
Simulation-wide correctness is verified through high-level tests.
*   **MMS (Method of Manufactured Solutions)**: Add more MMS cases to verify that solvers achieve the expected order of accuracy for the governing PDEs.
*   **Feature-Specific Inputs**: Create minimal `.ups` files that specifically target untested branches in the input file parsing and setup logic.
*   **Restart Consistency**: Ensure coverage for the restart logic by adding tests that compare a continuous run against one that is stopped and restarted.

## Phase 5: Automation and Maintenance
*   **CI Integration**: (Future) Integrate coverage reporting into the CI pipeline to ensure that new PRs do not decrease the coverage percentage.
*   **Coverage Targets**: Set incremental goals (e.g., "increase coverage by 5% every quarter").

## Appendix: Running Coverage Tests

To generate a baseline coverage report from the `dbg` directory, follow these steps:

1.  **Zero the counters**:
    Ensure you start with a clean slate by resetting all coverage counters.
    ```bash
    lcov --directory dbg --zerocounters
    ```

2.  **Execute Tests**:
    Run all unit tests to generate the raw coverage data (`.gcda` files).
    ```bash
    find dbg -name CTestTestfile.cmake -exec dirname {} \; | xargs -I {} ctest --test-dir {}
    ```

3.  **Generate .gcov files**:
    Create a `coverage` subdirectory and generate the `.gcov` files there.
    ```bash
    mkdir -p dbg/coverage
    cd dbg/coverage && find .. -name "*.gcda" -exec gcov {} +
    ```

4.  **Capture Coverage Data**:
    Collect the data from the build directory into a tracefile.
    ```bash
    cd dbg/coverage && lcov --directory .. --capture --output-file coverage.info --ignore-errors mismatch --preserve
    ```

5.  **Clean and Filter Results**:
    Filter out system headers, submodules, and the tests themselves to focus on the project source code.
    ```bash
    lcov --remove dbg/coverage/coverage.info '/usr/*' '*/submodules/*' '*/dbg/*' '*/UnitTests/*' '*/tests/*' --output-file dbg/coverage/coverage.cleaned.info --ignore-errors unused
    ```

5.  **Generate Report**:
    Generate a human-readable HTML report or a console summary.
    *   **HTML Report**:
        ```bash
        genhtml dbg/coverage/coverage.cleaned.info --output-directory dbg/coverage/html
        ```
    *   **Console Summary**:
        ```bash
        lcov --list dbg/coverage/coverage.cleaned.info
        ```
