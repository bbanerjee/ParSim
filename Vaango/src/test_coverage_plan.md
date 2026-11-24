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

## Clarifications for Unit Testing (Phase 2 & 3)
1.  **Scope**: Work through `src/Core` subdirectories alphabetically.
2.  **Depth**: Implement basic "Sanity/Construction" tests initially.
3.  **Dependencies**: Link against required internal libraries; use `gmock` only for high-complexity dependencies.
4.  **Integration**: Use the modern CMake pattern with `gtest_discover_tests`.
5.  **Tracking**: Work one subdirectory at a time, documenting progress in `src/progress_report.md`.

## Summary of Completed Unit Testing
### Core/Containers
- Established `UnitTests` directory and CMake integration.
- Added basic construction and sanity tests for `ConsecutiveRangeSet`, `Array1`, `RangeTree`, and `SuperBox`.
- Verified integration with `Vaango_Core_Math` for geometry dependencies in `SuperBox` tests.

### Core/DataArchive
- Established `UnitTests` directory and CMake integration.
- Implemented construction test for `DataArchive` with a mock `.uda` structure.
- Resolved issues with internal error exceptions during construction by providing a minimal `index.xml`.

### Core/Datatypes
- Established `UnitTests` directory and CMake integration.
- Added basic tests for `DenseMatrix` and `ColumnMatrix`, including construction, identity matrix generation, and matrix-vector multiplication.
- Verified polymorphic behavior through the `Matrix` base class interface.

### Core/Disclosure
- Established `UnitTests` directory and CMake integration.
- Added tests for `TypeDescription` lookup and `toString` functionality.
- Verified `TypeUtils` helper functions for fundamental types.
- Noted that types must be explicitly registered via `fun_getTypeDescription` before lookup in tests.

### Core/Exceptions
- Established `UnitTests` directory and CMake integration.
- Added basic tests for common exception types (`InternalError`, `ProblemSetupException`, `FileNotFound`, `ArrayIndexOutOfBounds`).
- Verified that error messages and exception types are correctly reported.

### Core/Geometry
- Established `UnitTests` directory and CMake integration.
- Added basic tests for `Point`, `Vector`, and `BBox`.
- Verified construction, arithmetic operations (Point + Vector, Point - Point), and bounding box extension/containment.

### Core/GeometryPiece
- Standardized `UnitTests` directory and updated CMake integration to use `gtest_discover_tests`.
- Added basic tests for `BoxGeometryPiece` and `SphereGeometryPiece`.
- Verified volume calculations, bounding box generation, and point-in-shape tests.

### Core/Grid
- Standardized `UnitTests` directory and CMake integration.
- Added basic construction tests for `Grid` and `Level`.
- Verified basic patch management within a level.

### Core/IO
- Established `UnitTests` directory and CMake integration.
- Added tests for `UintahZlibUtil` reading from compressed files.
- Verified token-based reading (`getString`, `getInt`, `getDouble`) and line reading with comment skipping.

### Core/Malloc
- Established `UnitTests` directory and CMake integration.
- Added basic tests for the custom `Allocator`.
- Verified `scinew`, `scinew[]`, and standard `malloc` functionality when using the project's allocator.

### Core/Math
- Established `UnitTests` directory and CMake integration.
- Added basic tests for `Matrix3` (construction, identity, determinant, inverse) and `MiscMath` (Min, Max, Clamp).
- Verified core numerical operations essential for geometry and physics calculations.

### Core/OS
- Standardized `UnitTests` directory and CMake integration.
- Added basic tests for `ProcessInfo` class.
- Verified retrieval of memory usage information and human-readable formatting.

### Core/Parallel
- Established `UnitTests` directory and CMake integration.
- Added basic tests for `Parallel` static accessors.
- Verified that `getRootProcessorGroup` throws `InternalError` when MPI is not initialized, ensuring robust error handling in serial environments.

### Core/ProblemSpec
- Established `UnitTests` directory and CMake integration.
- Added basic tests for `ProblemSpec` XML parsing from memory buffers.
- Verified attribute retrieval, block finding, and default value handling.

### Core/Util
- Established `UnitTests` directory and CMake integration.
- Added basic tests for `StringUtil` (conversion, case changing, splitting) and `Endian` (endianness detection, byte swapping).
- Verified core utility functions used throughout the codebase.

### Core/Lockfree (Skipped)
- Attempted to implement tests for `Lockfree::Pool`.
- Encountered persistent memory corruption errors attributed to interactions between the lockfree pool's memory management and the project's `sci-malloc` allocator.
- Decision made to skip this module for initial coverage efforts to maintain progress.

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
