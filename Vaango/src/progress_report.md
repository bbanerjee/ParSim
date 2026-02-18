# Vaango Unit Testing and Coverage Progress Report

## Summary
- **Current Coverage**: 7.1% (Lines)
- **Status**: Alphabetical unit testing pass through `src/Core` is COMPLETE.
- **Key Stability Fixes**: 
    - Disabled `Vaango_Core_Malloc` (sci-malloc) project-wide to resolve MPI initialization crashes and symbol lookup errors.
    - Modernized `src/Core/Malloc/new.cc` for C++20 `noexcept` compliance.
    - Fixed a critical memory corruption bug in `ArrayIndexOutOfBounds` exception.
    - Restored standard MPI/Environment initialization in unit tests.

## Completed Core Modules (100% Alphabetical Pass)
1.  **Containers**: `testArray1`, `testConsecutiveRangeSet`, `testRangeTree`, `testSuperBox`.
2.  **DataArchive**: `testDataArchive`.
3.  **Datatypes**: `testMatrices`.
4.  **Disclosure**: `testTypeDescription`.
5.  **Exceptions**: `testExceptions` (Fixed critical `char*` bug).
6.  **Geometry**: `testGeometry`.
7.  **GeometryPiece**: `testGeometryPieces`, `testLocalSDF`, `testDynamicSDFGeometry`.
8.  **Grid**: `testGrid`, `testMaterialManager`.
9.  **IO**: `testZlib`.
10. **Malloc**: `testAllocator` (Updated to handle `DISABLE_SCI_MALLOC`).
11. **Math**: `testMath`.
12. **OS**: `testDir`, `testProcessInfo`.
13. **Parallel**: `testParallel`, `testParallelInit`.
14. **ProblemSpec**: `testProblemSpec`.
15. **Util**: `testUtil`.

## Completed CCA Components
- **SimulationCommon**: `testSimulationCommon`, `testSimulationReductionVariable`.
- **DataArchiver**: `testDataArchiver`.
- **Schedulers**: `testSchedulerCommon`.
- **Regridder**: `testRegridderCommon`.
- **LoadBalancers**: `testLoadBalancerCommon`.
- **SolidReactionModel**: `testArrhenius`, `testModifiedArrhenius`, `testNthOrderModel`, `testPowerModel`.
- **ElasticModuliModels**: `testElasticModuli_Arena`, `testElasticModuli_MetalIso`, `testElasticModuli_NeuralNet`, `testElasticModuli_NeuralNet_Bulk`, `testElasticModuli_SupportVector`, `testElasticModuli_Tabular`, `testElasticModuli_Constant`, `testElasticModuli_Arenisca`, `testElasticModuli_ArenaMixture`, `testElasticModuli_Tabular_Bulk`, `testElasticModuli_Tabular_BulkPressure`.
- **EOSModels**: `testEOS_BorjaT`, `testPressure_Air`, `testPressure_Granite`, `testPressure_Water`, `testMieGruneisenEOS`, `testHyperElasticEOS`, `testDefaultHypoElasticEOS`, `testBorjaEOS`.
- **TabularModels**: `testTabularData`, `testTabularEOS`, `testTabularPlasticity`, `testTabularPlasticityCap`, `testNeuralNetTabularPlasticity`.

## Blockers / Skipped
- **Core/Lockfree**: Skipped due to persistent allocator conflicts with `std::allocator`.

## Next Steps
1.  Extend testing to remaining `CCA` components.
2.  Maintain `DISABLE_SCI_MALLOC` for all builds to ensure CI stability.
