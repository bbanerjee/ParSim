# Design Document: Hybrid DEM-MPM Coupling in Vaango

## Overview
This document outlines the design for integrating Discrete Element Method (DEM) capabilities into the existing MPM framework within Vaango. The primary goal is to support simulations of jointed rock masses and other systems where continuum behavior (MPM) and discrete, arbitrary-shaped interactions (DEM) coexist.

---

## 1. Hybrid DEM-MPM Approach (Particle-Based)

In this mode, DEM particles are treated as special MPM particles that bypass continuum stress-divergence calculations and instead use pair-wise contact laws.

### Data Structures
- **Rotational State**: New `ParticleVariable`s for `p.angularVelocity`, `p.torque`, and `p.orientation`.
- **Geometric Properties**: `p.radius` for spherical contact.
- **Material Identity**: `MPMMaterial` extension to include an `is_dem_material` flag.

### Task Graph Workflow
1. **Neighbor Finding**: Utilize the MPM background grid for $O(N)$ binning. Particles in neighboring cells (within $2 \times R_{max}$) are identified as potential contact pairs.
2. **ComputeDEMForces Task**: 
   - Iterates through particles in cells.
   - Performs collision detection: $\|x_i - x_j\| < (r_i + r_j)$.
   - Applies contact laws (e.g., Hertzian, Linear Spring-Dashpot).
   - Accumulates resulting forces into `p.externalForce` and torques into `p.torque`.
3. **Momentum Transfer**: `p.externalForce` is naturally transferred to grid nodes during standard MPM interpolation, implicitly coupling the discrete and continuum phases.
4. **Rotational Integration**: A new task updates angular kinematics:
   - $\omega^{n+1} = \omega^n + (\text{Torque} / \text{Inertia}) \Delta t$

---

## 2. Arbitrary Shaped Objects (SDF-Based)

For non-discretized, irregular rock blocks, we use Signed Distance Fields (SDF) to model geometry and detect contacts.

### Representation: `DynamicSDFGeometry`
- Inherits from `GeometryPiece`.
- Holds the rigid body state (Position, Orientation, Linear/Angular Velocity).
- Stores a local SDF grid representing the object surface ($\phi < 0$ inside, $\phi > 0$ outside).

### Contact Detection Mechanism
- **DEM-to-DEM**: Contact occurs if points on the surface of Block B have a negative value in the SDF of Block A. The penalty force acts along $\nabla \phi_A$.
- **MPM-to-DEM**: MPM particles (treated as points) check their position against the `DynamicSDFGeometry` SDF. If inside, a penalty force is applied to the particle and an equal/opposite force/torque to the rigid block.

### Advantages for Jointed Rock Masses
- **Asperity Resolution**: Captures joint roughness via the SDF without requiring millions of surface particles.
- **Efficiency**: Large rock blocks are handled as single rigid bodies (or coarse meshes), significantly reducing the computational load compared to full discretization.

---

## 3. Implementation Plan

### Phase 1: Core Definitions
- Update `MPMLabel.cc/h` with new particle labels (`p.angularVelocity`, `p.radius`, etc.).
- Update `MPMFlags.cc/h` and `MPMMaterial.cc/h` to parse DEM-specific parameters from the `.ups` file.

### Phase 2: Hybrid Task Integration
- Implement `scheduleComputeDEMForces` in `SerialMPM.cc`.
- Implement `scheduleIntegrateDEMRotation` to handle angular momentum.

### Phase 3: Geometry Piece Promotion
- Create `Core/GeometryPiece/DynamicSDFGeometry`.
- Implement SDF generation from `TriGeometryPiece` (STL/OBJ input) using a Fast Marching Method.

### Phase 4: Verification
- Validate with simple sphere-sphere collisions (analytical solution).
- Validate DEM block-to-MPM bed interaction (e.g., a rock falling into soil).

---

## 4. Implementation Summary (Completed)

### Key Implementation Details:

1. **SDF Infrastructure**:
    - **`Core/GeometryPiece/LocalSDF`**: A new utility class that performs trilinear interpolation on a local Cartesian grid to provide distance and gradient (normal) information. It also handles world-to-local coordinate transformations.
    - **`Core/GeometryPiece/DynamicSDFGeometry`**: A new geometry piece that wraps a `TriGeometryPiece` (STL/OBJ). It pre-calculates an SDF upon initialization, enabling complex shapes to be used as rigid bodies.

2. **Hybrid Task Integration in `SerialMPM`**:
    - **`scheduleComputeDEMForces`**: This task is inserted early in the time-step. It uses ghost-cell data to detect contacts between particles.
    - **Dual Contact Logic**:
        - **Spherical DEM**: Uses a Linear Spring-Dashpot model for pair-wise contact between particles marked as `is_dem_material`.
        - **SDF Rigid Bodies**: Particles representing `DynamicSDFGeometry` pieces act as rigid blocks. The code queries the SDF to detect overlap with continuum MPM particles or other blocks, applying penalty forces and torques.
    - **`scheduleIntegrateDEMRotation`**: A new integration task that updates angular velocity and orientation based on computed torques and the body's inertia tensor.

3. **Data Management and Relocation**:
    - **State Variables**: Added `p.angularVelocity`, `p.torque`, `p.orientation`, `p.radius`, and `p.inertiaTensor` to the permanent particle state.
    - **Automatic Parallelism**: These blocks are represented as single particles in Vaango’s `DataWarehouse`, meaning they automatically relocate across MPI ranks and patches as they move through the domain.
    - **Efficient Creation**: The `ParticleCreator` was modified to recognize `DynamicSDFGeometry` and create exactly one rigid-body particle at its center, rather than discretizing its interior.

4. **User Control**:
    - The features are activated via the `<enable_dem>true</enable_dem>` flag in the `MPM` section of the input file.
    - Individual materials can be designated as discrete and/or SDF-based using `<is_dem_material>true</is_dem_material>` and `<is_sdf_based>true</is_sdf_based>`.

### Summary of Added/Modified Files:
- **New**: `Core/GeometryPiece/LocalSDF.h/.cc`, `Core/GeometryPiece/DynamicSDFGeometry.h/.cc`, `dem_mpm_coupling.md`
- **Modified**: `MPMLabel.h/.cc`, `MPMFlags.h/.cc`, `MPMMaterial.h/.cc`, `SerialMPM.h/.cc`, `ParticleCreator.h/.cc`, `GeometryPieceFactory.cc`, `Core/GeometryPiece/CMakeLists.txt`

---

## 5. Symmetrical Multi-Material Task Setup

To ensure momentum conservation (Newton's Third Law) and improve efficiency, the DEM force computation is implemented as a MaterialSet-wide task.

### Key Refinements:
1. **MaterialSet-wide Scheduling**: The `computeDEMForces` task is scheduled once for the entire `MaterialSet` of MPM materials. It uses `modifies` access for `p.externalForce` and `p.torque` across all materials.
2. **Newton's Third Law**: Forces are calculated for a pair $(i, j)$ once. $\mathbf{F}$ is applied to particle $i$ and $-\mathbf{F}$ is applied to particle $j$.
3. **Ghost Particle Handling**: To maintain consistency in parallel simulations:
    - Particle $i$ (the "owner" in the current loop) always receives the force.
    - Particle $j$ receives the equal and opposite force **only if it is a real particle** on the current patch.
    - If particle $j$ is a ghost, the corresponding force will be calculated and applied by the patch that "owns" it.
4. **Performance**: Reduces collision checks by roughly 50% by avoiding redundant pair calculations.

---

## 6. UPS_SPEC Updates

The following tags have been added to the input specification files in `StandAlone/inputs/UPS_SPEC` to enable validation and documentation of the new hybrid features.

### Global Flags (`mpm_spec.xml`)
- `<enable_dem>`: (Boolean) Enables the Discrete Element Method infrastructure.

### Material Flags (`mpm_spec.xml`)
- `<is_dem_material>`: (Boolean) Designates a material as a discrete phase.
- `<is_sdf_based>`: (Boolean) Specifies that the material uses SDF-based rigid body contact.
- `<radius>`: (Double) The characteristic radius of the discrete particles.
- `<kn>`, `<kt>`: (Double) Normal and tangential contact stiffness.
- `<mu>`: (Double) Friction coefficient.
- `<gamma>`: (Double) Viscous damping coefficient.

### Geometry Pieces (`ups_spec.xml`)
- `<dynamic_sdf>`: A new geometry piece type.
    - Requires `<res>`: (Vector) Grid resolution for the SDF.
    - Requires `<tri>`: (Block) The underlying surface mesh definition.
    - Can be used inside `<union>`, `<difference>`, and `<intersection>`.

---

## 7. Handling Large DEM Objects Across MPI Boundaries

A potential issue with the "One-Particle-per-Body" approach is that contact between large DEM blocks in adjacent MPI patches may be missed if the particle (representing the centroid) is outside the ghost-cell range of the neighbor patch.

### Suggested Solutions:

1. **Radius-Dependent Ghosting Range**:
   - In `scheduleComputeDEMForces`, dynamically calculate the required ghost cell depth based on the maximum radius ($R_{max}$) of all discrete/SDF materials: `num_ghost = ceil(R_max / dx)`.
   - *Pros*: Simple to implement using existing infrastructure.
   - *Cons*: High communication overhead if a few objects are very large compared to the grid size.

2. **Representative "Slave" Particles**:
   - Instead of 1 particle at the center, represent each large block with a sparse "cloud" of particles (e.g., center + 8 bounding box corners).
   - These slave particles carry the same rigid body ID and trigger Vaango's standard ghosting/relocation logic in every patch the block physically touches.
   - *Pros*: Reuses optimized ghosting paths; communication is proportional to the object's surface area/span rather than a global $R_{max}$.

3. **Broad-Phase Patch-Object Registry**:
   - Implement a lightweight global broad-phase task that maps each large `DynamicSDFObject` to the patches it overlaps.
   - Explicitly request the object's state particle as a requirement for all overlapping patches, regardless of centroid distance.
   - *Pros*: Most efficient for memory and communication.
   - *Cons*: Requires custom modifications to the task-dependency logic in the `Scheduler`.

### Recommendation for Rock Masses:
For jointed rock masses where blocks are large and irregular, **Solution 2 (Slave Particles)** is recommended as it provides the best balance between implementation complexity and parallel efficiency.

---

## 8. Implementation of Slave Particles (Geometric Proxies)

The "Slave Particles" approach has been implemented to ensure that large rigid bodies are correctly detected for contact even when their centroids are far from patch boundaries.

### Key Implementation Details:
1. **Distributed Representation**: For each `DynamicSDFGeometry`, the `ParticleCreator` now generates **9 particles**:
    - **1 Master**: Located at the centroid, carries the full mass and inertia of the block.
    - **8 Slaves**: Located at the corners of the bounding box, acting as "geometric proxies" with zero mass.
2. **Rigid Body Identification**:
    - A new `p.rigidBodyID` label was added to the permanent particle state.
    - All 9 particles belonging to the same block are assigned the same unique `rbID` (derived from material index and geometry object index).
3. **Contact Filtering**:
    - `SerialMPM::computeDEMForces` was updated to pre-fetch `p.rigidBodyID`.
    - A collision check `if (rbID_i == rbID_j) continue;` ensures that master and slave particles belonging to the same rigid body do not interact with each other.
4. **MPI Visibility**:
    - Because the 8 slave particles are at the object's extremities, they trigger Vaango's standard ghost-cell and relocation logic in any patch the block physically touches.
    - This ensures that neighbor patches "see" the rigid block and can calculate contact forces between their local particles and the block's proxy particles.
5. **Mass Handling**:
    - Slave particles are assigned `p.mass = 0.0`. This distinguishes them from the master particle during dynamics integration and ensures that they do not contribute extra mass to the system or the MPM grid.

---

## 9. Unit Testing

A comprehensive unit test for the hybrid DEM-MPM implementation has been added to the project.

### Test Details:
- **Location**: `CCA/Components/MPM/UnitTests/testHybridDEM.cc`
- **Test Case: `CollisionTest`**:
    - Programmatically generates a `.ups` file using `libxml2`.
    - Configures a simulation with a discrete material (`is_dem_material = true`).
    - Places two spherical particles such that they overlap initially.
    - Initializes the full Vaango simulation environment (`SimulationController`, `Scheduler`, etc.).
    - Runs the simulation for 2 timesteps.
    - This verifies the entire integration chain: from XML parsing of new tags, through particle creation (including Master/Slave logic), to the execution of the symmetrical `computeDEMForces` task.
- **Build System**: The `UnitTests` directory is integrated into the `CCA/Components/MPM` build via `add_subdirectory`.

### Core Geometry Unit Tests:
- **Location**: `Core/GeometryPiece/UnitTests/`
- **`testLocalSDF.cc`**:
    - Verifies coordinate transformations between world and local space.
    - Verifies trilinear interpolation of distance values and gradients.
- **`testDynamicSDFGeometry.cc`**:
    - Verifies object creation from XML specification.
    - Checks exception handling for missing mesh files.

---

## 10. Alternative: Scheduler-Based Large Object Handling

Modifying the Scheduler to handle large DEM blocks involves moving away from fixed-range spatial ghosting towards dynamic, object-aware dependency tracking.

### Architectural Breakdown:

1. **Broad-Phase Patch-Object Map**:
   - New Task: `scheduleMapBlocksToPatches` (executed at the start of every timestep).
   - Action: Calculates intersection of each block's oriented bounding box with grid patches.
   - Output: A global registry mapping `BlockID -> List<PatchID>`.

2. **Custom Dependency Injection**:
   - Mechanism: Override grid-based intersection logic. If Patch A contains Block B, and the map says Block B overlaps Patch C, create an explicit dependency: *Patch C requires Particle B from Patch A*.
   - Code Change: Add `task->needsParticle(label, particleID)` to the `Task` class.

3. **Modifying OnDemandDataWarehouse**:
   - Current State: Packs particles within spatial ghost range.
   - Modified State: Update to pack "Explicitly Requested Particles" (e.g., Block B) for Patch C, even if the centroid is outside the ghost range.

4. **Updating the Relocator**:
   - Problem: Ensure correct "Master" handoff and notification of overlapping patches when a block's center moves between patches.

### Comparison:
- **Slave Particles**: Recommended for jointed rock masses (thousands of blocks). Reuses existing parallel machinery.
- **Scheduler Modification**: Best for simulations with a few massive objects (e.g., tectonic plates) where slave particle overhead is prohibitive.

---

## 11. Detailed Architecture: Hybrid Structured+Graph Scheduler

To handle large objects efficiently without sacrificing the performance of regular grids, the Scheduler can be evolved into a hybrid system.

### Core Concept: Dual Dependency Engines
The Scheduler will operate using two distinct dependency engines that function in parallel:
1.  **Structured Engine (Mesh Data)**: Handles 99% of the simulation data (continuum particles, grid nodes). Uses the existing $O(1)$ `Task::requires(Ghost::AroundNodes)` logic, which infers dependencies implicitly from grid indices.
2.  **Graph Engine (Large Objects)**: Handles the remaining 1% (large rigid bodies). Uses explicit $O(E)$ graph edges created dynamically based on object-patch overlaps.

### Detailed Technical Changes:

#### 1. Task Dependency Representation
- **Current**: `Task::Dependency` holds `Patch*`, `VarLabel*`, `GhostType`.
- **New**: Extend `Task::Dependency` to include `std::set<long64> particleIDs`.
- **Flagging**: Introduce a dependency type flag: `Dependency::Type = Global | Spatial | ExplicitParticle`.

#### 2. DetailedTasks Graph Construction (The "Compile" Phase)
- **Standard Phase**: Iterate through patches to find spatial neighbors using the existing efficient logic.
- **Overlay Phase**: Iterate through `DynamicSDFObjects`.
    - For Object $O$, calculate its bounding box $B$.
    - Identify all patches $P_{overlap}$ that intersect $B$.
    - For every Task $T$ that consumes Object $O$:
        - Create explicit dependency edges from the *Current Owner Patch* of $O$ to all patches in $P_{overlap}$.
        - Mark these edges as "Sparse Data Transfer".

#### 3. DataWarehouse (The "Runtime" Phase)
- **Packing**: Modify `OnDemandDataWarehouse::sendMPI`.
    - Add a "Sparse Particle Buffer" to the packing logic.
- **Execution**: When processing a dependency edge:
    - If **Spatial**: Pack the grid region (existing logic).
    - If **Explicit**: Perform a serialized binary copy of *only* the specific particle indices associated with the requested ID.

#### 4. Load Balancing (Dynamic)
- **Challenge**: Large objects introduce "Long-Range Dependencies," increasing communication latency.
- **Cost Model Update**: Update the load balancer to add a penalty weight to patches hosting the Master Particles of large objects. This accounts for the extra serialization and communication overhead required to serve the object to multiple remote patches.

### Efficiency Strategy
- **Bypass**: The graph engine is *only* engaged for particles flagged as `is_large_object`. Standard grid variables (`g.mass`, `g.velocity`) completely bypass the graph checks.
- **Caching**: The `Patch-to-Object` map is cached and only recalculated during `Relocation` intervals, preventing overhead during every micro-step.

---

## 12. Stability Conditions for Hybrid Integration

In explicit time integration, numerical stability depends on resolving the highest frequency present in the system. While standard MPM is governed by the CFL condition (wave speed on the grid), DEM is governed by the contact stiffness.

### DEM Stability Limit:
For a linear spring-mass system (as used in our `computeDEMForces` task), the critical time step is:
$$\Delta t_{crit} = \pi \sqrt{\frac{m}{k_n}}$$

### Implementation:
1. **Task Update**: `SerialMPM::actuallyComputeStableTimestep` has been updated to iterate through all discrete particles.
2. **Calculation**: For each master particle (where $m > 0$), it calculates:
   $$\Delta t_{particle} = \beta \sqrt{\frac{m}{k_n}}$$
   where $\beta = 0.2$ is a safety factor to account for tangential forces and multi-contact scenarios.
3. **Coupling**: The final simulation time step is the minimum of:
   - Grid-based MPM CFL limit (from constitutive models).
   - Particle-based DEM stability limit (from this task).
   - Any user-defined limits in the `.ups` file.

This ensures that the simulation remains stable even when using high contact stiffness values typical of jointed rock masses.

---

## 13. Conditional Registration and Allocation

To maintain compatibility with standard MPM simulations and avoid task graph compilation errors during relocation, all DEM-specific logic has been made conditional.

### Mechanism:
1. **Registration**: In `ParticleCreator::registerPermanentParticleState`, labels like `p.rigidBodyID`, `p.angularVelocity`, etc., are only added to the `particle_state` vectors if `d_flags->d_enableDEM` is true.
2. **Allocation**: Memory for these variables is only allocated in `ParticleCreator::allocateVariables` when DEM is enabled.
3. **Initialization**: The ` एक्चुअलीInitialize` and `refine` tasks in `SerialMPM` only attempt to populate these variables if `enable_dem` is active.

### Benefits:
- **Zero Overhead**: Standard MPM simulations do not allocate memory for or communicate DEM-specific data.
- **Robustness**: Fixes the "Failed to find comp for dep" error in the `Scheduler` by ensuring that only variables being produced by active tasks are requested for relocation.

---

## 14. Robust Discrete Object Handling and Stability

Recent refinements ensure that the Slave Particles approach is robust against Vaango's internal stability mechanisms and numerical edge cases.

### 1. Particle Deletion Protection
Vaango typically deletes particles whose mass falls below a threshold (`d_minPartMass`) or whose temperature becomes unphysical. Since Slave particles are designed to have zero mass, they would normally be deleted.
- **Fix**: Modified `SerialMPM.cc` to exempt materials flagged with `isDEMMaterial()` from mass-based deletion.
- **Benefit**: Ensures that Slave particles persist as geometric proxies throughout the simulation, regardless of their mass.

### 2. Isolated Grid Contributions
To prevent zero-mass particles from interfering with the continuum phase calculations:
- **Change**: Updated `interpolateParticlesToGrid` to skip any particle with $m \le 0$.
- **Result**: Slave particles do not contribute to grid mass, volume, or momentum, while the Master particle (carrying the total object mass) provides the correct physical coupling.

### 3. Singular Inertia Handling
Numerical errors can occur if rotation integration is attempted on particles with zero mass or radius (resulting in a singular inertia matrix).
- **Fix**: Added a determinant check in `integrateDEMRotation`. Integration is only performed if $\det(\mathbf{I}) > 0$.
- **Fallback**: Slave particles skip rotation integration (as they are geometrically tied to the Master), preventing $NaN$ propagation.

### 4. Geometry-Aware Initialization
- **Cylinders/Spheres**: `ParticleCreator` now inherits the actual radius from the `GeometryPiece` rather than relying on global material defaults.
- **Boxes**: A representative radius is calculated as $0.5 \times \text{smallestSide}$.
- **Ordering**: Mass, volume, and radius are fully initialized *before* the inertia tensor is calculated, ensuring consistency.

---

## 15. Distributed Contact Model and Refined Collision Detection (February 2026)

To resolve issues where large rigid bodies (DEM) only interacted with continuum particles at a single point, the contact logic and supporting math libraries were significantly updated.

### 1. Distributed Contact Logic in `SerialMPM`
The original "centroid-to-particle" contact model was replaced with a **Distributed Contact Model**:
- **Mechanism**: For every continuum particle $i$, the code now identifies the **closest** discrete surface particle $j$ belonging to a rigid body.
- **Surface Reaction Forces**: Reaction forces are applied directly to the surface particle $j$ rather than being concentrated at the rigid body's master particle (centroid).
- **Physical Consistency**: While forces are applied at the surface to ensure correct momentum transfer at the grid level, they are also aggregated (along with torques) at the **Master particle** to drive the rigid body's global translation and rotation.
- **Visualization**: This change ensures that UDA outputs for `p.externalforce` correctly show a distributed interaction across the entire contact surface.

### 2. Direction-Dependent Effective Radius
Continuum particles in MPM are typically treated as points in DEM contexts. To ensure contact is detected the moment a particle's volume overlaps with a DEM surface:
- **Calculation**: An effective radius $R_{\text{eff}}$ is calculated based on the particle's size vectors ($\mathbf{V}_1, \mathbf{V}_2, \mathbf{V}_3$) from `p.size` and the contact normal $\mathbf{n}$:
  $$R_{\text{eff}} = 0.5 \times (|\mathbf{V}_1 \cdot \mathbf{n}| + |\mathbf{V}_2 \cdot \mathbf{n}| + |\mathbf{V}_3 \cdot \mathbf{n}|)$$
- **Impact**: This allows for accurate contact detection for non-spherical grid-aligned particles, resolving gaps between the "point" centers and the actual continuum surface.

### 3. Math Library Extensions (`Matrix3`)
To support the $R_{\text{eff}}$ calculation, the `Matrix3` class was extended:
- **New Method**: `getColumn(int j)` was implemented as an inline method in `Matrix3.h`.
- **Functionality**: Returns a `Vector` representing the $j$-th column of the matrix, enabling efficient access to the particle's basis vectors stored in the `p.size` tensor.

### 4. Observability and Debugging
- **`DEM` DebugStream**: A new `DebugStream` named `DEM` was added to `SerialMPM.cc`.
- **Activation**: Users can enable detailed per-particle collision reports (including `phi`, `overlap`, and `totalForce`) by setting the environment variable `SCI_DEBUG="DEM:+"`.
- **Task Requirements**: `scheduleComputeDEMForces` was updated to explicitly require `p.size` to support the new radius logic.

---

## 16. Surface Dummy Particles for Visualization (February 2026)

To provide a clearer visual representation of discrete objects during post-processing, a system for generating and synchronizing surface "dummy" particles was implemented.

### 1. Geometry-Based Particle Generation
The `ParticleCreator` was updated to generate a cloud of particles on the surface of discrete objects:
- **Cylinders**: A mesh of particles is distributed along the side surface and on both end caps.
- **Spheres**: A longitude/latitude grid of particles is generated on the surface.
- **Properties**: All surface dummy particles are initialized with **zero mass**, **zero volume**, and **zero radius**. This ensures they remain purely visual and do not contribute to the system's mass, energy, or grid-level momentum.

### 2. Kinematic Synchronization
To ensure these particles accurately represent the rigid body's motion, they are synchronized with the master particle (centroid) at every timestep:
- **Position Update**: Dummy particles are translated based on the master's position update.
- **Rotational Mapping**: If the rigid body has angular velocity, the position of each dummy particle is updated using the rigid body rotation formula: $\mathbf{v}_p = \mathbf{v}_{master} + \boldsymbol{\omega} \times \mathbf{r}_{local}$.
- **State Inheritance**: Orientation and angular velocity are carried forward from the master particle to all associated dummy particles in the `integrateDEMRotation` task.

### 3. Contact Filtering
To maintain physical accuracy, the `computeDEMForces` task was modified to explicitly ignore particles with $m \le 0$. This ensures that the dummy particles do not participate in contact detection or generate spurious reaction forces, while the distributed contact model correctly handles interactions via the legitimate surface-proxy particles.

---

## 18. Distinct Treatment of Rigid Materials (February 2026)

To align with the physical requirements of rigid-body MPM, the implementation was updated to distinguish between `is_dem_material` and `is_rigid` materials.

### 1. Rigid Material Discretization
- **Logic**: Materials flagged as `is_rigid` are discretized using the standard grid-based approach (multiple particles per cell).
- **Master Reference**: A single master particle is maintained at the centroid to hold global kinematics, but all particles carry mass.

### 2. Grid-Based Contact
- **MPM Rules**: Rigid materials are excluded from pairwise DEM contact checks. They interact via standard MPM grid-based contact (specified, penalty, or friction contact).
- **Rigid Enforcement**: Particle velocities are updated from the grid normally and then "rigidified" by projecting them onto the translation/rotation of the body's centroid.

### 3. Discrete vs. Rigid Summary
| Feature | Discrete Material (`is_dem_material`) | Rigid Material (`is_rigid`) |
| :--- | :--- | :--- |
| **Discretization** | 1 Master + N Surface Slaves | Full MPM Grid Discretization |
| **Contact Model** | Pairwise DEM (Spring-Dashpot) | Standard MPM (Grid-Based) |
| **Mass Handling** | Only Master carries mass | All particles carry mass |
| **Internal Force** | N/A | Zero (handled by CM) |
| **Kinematics** | Driven by pairwise forces | Driven by grid-contact forces |

---

## 19. Integration of Rigid-Standard MPM Contact (February 2026)

To correctly handle contact between rigid MPM tools and standard MPM materials, the force extraction and kinematic enforcement were refined.

### 1. Contact Force Extraction
- **Logic**: In `interpolateToParticlesAndUpdate`, the total external force (including contact) is now explicitly calculated for all particles:
  $$\mathbf{f}_{ext}^{total} = m \frac{\mathbf{v}_{new} - \mathbf{v}_{old}}{\Delta t} - \mathbf{f}_{body}$$
- **Visualization**: This value is stored in `p.externalforce`, ensuring that UDA outputs show the correct interaction forces during both DEM and MPM-style contacts.
- **Consistency**: By populating `p.externalforce` with the total resultant force, the subsequent timestep's `interpolateParticlesToGrid` naturally provides the correct `g.externalforce` for rigid body acceleration.

### 2. Resultant Force Acceleration for Rigid Bodies
- **Global Motion**: In `computeAndIntegrateAcceleration`, materials flagged as `is_rigid` now compute a single, uniform acceleration based on the sum of all nodal external, body, and internal forces.
- **Kinematic Constraints**: The "rigidification" pass in `interpolateToParticlesAndUpdate` ensures that all particles in a rigid body follow the master centroid's translation and rotation, even when interacting via grid-based contact rules.

### 3. Separation of Concerns
- **DEM (`is_dem_material`)**: Continues to use pairwise contact laws.
- **Rigid MPM (`is_rigid`)**: Uses standard MPM grid contact but enforces rigid kinematics and aggregates resultant forces.

---

## 20. Corrected Physical Representation of Rigid Materials (February 2026)

Following a review of the physical requirements for rigid bodies in hybrid simulations, the implementation was corrected to distinguish between point-proxy discrete bodies and grid-discretized rigid bodies.

### 1. Unified Discretization
- **Requirement**: Rigid materials must be discretized identically to deformable materials to ensure standard MPM contact rules are valid.
- **Change**: Reverted the master/slave sparse discretization for `is_rigid` materials. They now generate a full cloud of particles.
- **Mass Distribution**: Every particle in a rigid material carries mass, providing accurate momentum and volume contributions to the grid.

### 2. Grid-Driven Resultant Kinematics
- **Mechanism**: The rigid body's motion is now driven by the **resultant** of all grid-level forces (external, body, and contact).
- **Consolidation**: In `computeAndIntegrateAcceleration`, a single uniform acceleration $\mathbf{a}_{rigid} = \sum \mathbf{F} / \sum M$ is calculated from all nodes associated with the rigid material.
- **Enforcement**: This acceleration is used to update the master reference particle, and a subsequent pass in `interpolateToParticlesAndUpdate` synchronizes all associated particles to this rigid motion.

### 3. Verification of Contact Force Output
- **Explicit Population**: By extracting the contact force as the difference between the grid-interpolated momentum change and prescribed body forces, the `p.externalforce` variable now correctly visualizes the mechanical interaction between rigid tools and deformable media.
- **Standard Contact Compatibility**: This approach ensures that standard contact models (Friction, Specified Body, Penalty) work out-of-the-box with rigid materials while maintaining non-deformability.


---

## 21. Unphysical Rotation of DEM Materials (February 2026)

Based on an examination of the input file and the source code, the rotation of the DEM piston is likely caused by an extremely underestimated moment of inertia.

###  Analysis:
- **Incorrect Inertia Calculation**: In ParticleCreator.cc, the moment of inertia for discrete materials is hardcoded using the formula for a solid sphere: $I = 0.4 \times \text{mass} \times \text{radius}^2$.  For a BoxGeometryPiece, the code sets the "radius" to be half of the smallest side of the box.  In your simulation, the piston box dimensions are approximately $1.0 \times 0.0125 \times 0.1$. The "radius" is thus set to $0.00625$.  The resulting rotational inertia ($I \approx 1.56 \times 10^{-5} \times m$) is used for all axes, whereas the actual inertia for a rectangle of width 1.0 should be.
- **Symmetry Breaking**: In any numerical simulation, there are small asymmetries due to floating-point precision and the discrete nature of the MPM particles in the deformable HMX cylinder. When the piston descends and contacts the cylinder, these small asymmetries produce a net torque.
- **Amplification**: Because the rotational inertia is so unnaturally small, even a negligible torque results in a massive angular acceleration, causing the piston to spin wildly during the simulation.

###  Recommendation:
  The moment of inertia calculation in ParticleCreator::initializeParticle should be updated to account for the actual geometry type and dimensions of the GeometryPiece.  Specifically for a BoxGeometryPiece, it should use the standard formula for a rectangular cuboid: $I_{yy} = \frac{1}{12} m (a^2 + c^2)$, etc., where $a, b, c$ are the side lengths.

---

## 22. Geometry-Aware Rotational Inertia and Dynamics (February 2026)

To resolve unphysical rotations caused by underestimated inertia, the `GeometryPiece` hierarchy was extended to provide exact physical properties, and the rotational integration in `SerialMPM` was refined.

### 1. Physical Property Extensions in `GeometryPiece`
The base `GeometryPiece` class was updated to include virtual methods for calculating exact volume, center of mass, and the moment of inertia tensor.
- **Exact Formulas**: Implemented geometry-specific inertia tensors for all core shapes:
    - **`BoxGeometryPiece`**: Rectangular cuboid formulas.
    - **`SphereGeometryPiece`**: Solid sphere formulas.
    - **`CylinderGeometryPiece` / `ConeGeometryPiece`**: Formulas for oriented cylinders and frustums.
    - **`TriGeometryPiece`**: Summation of tetrahedron contributions for closed triangular meshes.
    - **`Union` / `Difference`**: Calculated using the **Parallel Axis Theorem** to correctly shift and combine child inertia tensors.
- **Orientation Support**: Inertia tensors are calculated in the object's local frame and then rotated into the global basis using the object's initial orientation.

### 2. Accurate Particle Initialization
The `ParticleCreator::initializeParticle` method was refactored to use these new geometric queries:
- **Mass and Volume**: Now uses the exact `piece->volume()` and `piece->getCenter()`.
- **Inertia Tensor**: The initial `p.inertiaTensor` is populated by querying the `GeometryPiece` and scaling by the material's initial density. This ensures that thin or elongated objects (like the piston) have the correct resistance to rotation about all axes.

### 3. Inertia Tensor Evolution
In `SerialMPM::integrateDEMRotation`, the moment of inertia tensor is now evolved as the rigid body rotates:
- **Mechanism**: The rotation increment $\Delta \mathbf{R} \approx \mathbb{I} + \boldsymbol{\Omega} \Delta t$ is used to update the tensor:
  $$\mathbf{I}_{new} = \Delta \mathbf{R} \cdot \mathbf{I}_{old} \cdot \Delta \mathbf{R}^T$$
- **Consistency**: This ensures that as an object rotates, its resistance to further rotation correctly follows its spatial orientation.
- **Synchronization**: The updated inertia tensor and orientation are correctly propagated from the Master particle to all associated Slave particles (geometric proxies).

### 4. Box Surface Visualization
To support better visual inspection of box-shaped rigid bodies, a new `createBoxSurfacePoints` helper was added to `ParticleGeometryHelpers`. This generates a grid of visual-only particles on all six faces of a `BoxGeometryPiece`, providing a clear boundary for contact verification in post-processing.

## 23. Comparison with State-of-the-Art (Liu et al., 2018)

A comparison with the coupling algorithm presented by Liu et al. (2018) highlights the unique design choices in Vaango's hybrid DEM-MPM implementation.

### 1. Representation of DEM Bodies
*   **Liu et al.**: Each DEM block is represented by **nine "visual material points"** (vertices, midpoints, and center). These points are used to project mass and momentum onto the MPM background grid.
*   **Vaango**: DEM bodies are represented by a **Master-Slave particle system**. Instead of a fixed set of nine points, Vaango utilizes arbitrary **Geometry Pieces** (Box, Sphere, Cylinder, TriMesh) and can support an arbitrary number of slave (proxy) particles to define the surface. Furthermore, Vaango uses a **Signed Distance Field (SDF)** to represent the body's boundary precisely.

### 2. Contact Detection Mechanism
*   **Liu et al.**: The interaction between the granular flow (MPM) and the blocks (DEM) is handled via the **standard MPM grid-based multi-material contact algorithm**. Collisions are detected only if material points from both bodies contribute to the same grid node.
*   **Vaango**: Interaction is handled through a **direct Particle-SDF interaction** in the `computeDEMForces` task. For every standard MPM particle, we perform an SDF check against the DEM body's geometry. This allows for surface-to-surface contact that is independent of the grid resolution for the collision detection phase itself.

### 3. Force Calculation and Projection
*   **Liu et al.**: The contact force is derived from the **acceleration of grid nodes** and then projected back to the DEM vertices as a body force.
*   **Vaango**: Contact forces (Normal and Tangential) are calculated using a **penalty method** (stiffness $k_n, k_t$ and damping $\gamma$) directly at the contact point on the SDF surface. These forces and the resulting torques are applied directly to the **Master particle** of the DEM body, which governs its global rigid-body dynamics.

### 4. Integration of Rotational Motion
*   **Liu et al.**: Focuses on "deformable" blocks where the motion is governed by the displacement of vertices and an internal stiffness matrix.
*   **Vaango**: Implements explicit **rigid-body rotational dynamics**. It tracks and integrates the orientation (rotation matrix), angular velocity, and the orientation-dependent inertia tensor for each DEM body.

### 5. Distributed Contact (Recent Fix)
*   The recent correction in `SerialMPM::computeDEMForces` (now in `DEMTasks::computeDEMForces`) ensures that **every individual MPM particle** is checked against the DEM body. This avoids the "single-point contact" approximation that often occurs in grid-based methods (like Liu's) where multiple particles might be collapsed into a single nodal momentum contribution before contact is resolved.

## 24. Comparison with State-of-the-Art (Ren et al., 2022)
Based on an examination of the paper by Ren et al. (2022) and the current
implementation in Vaango, the primary differences in the coupling algorithm
are as follows:

### 1. Representation of Discrete Bodies
* Ren et al.: Uses the spheropolygon (sphero-polyhedra in 3D) technique to
    represent non-spherical particles. This involves "eroding" a polygon by a
    radius $r_{sphero}$ and then "dilating" it, resulting in a shape with
    rounded corners.
* Vaango: Uses Signed Distance Fields (SDF) combined with a Master-Slave
    particle system. Vaango can represent arbitrary shapes (Box, Cylinder,
    Sphere, TriMesh) using an analytical or grid-based SDF. The surface is
    explicitly populated with visual/proxy particles for interaction and
    visualization.

### 2. Multi-Phase Framework
* Ren et al.: Implements a two-phase (solid and liquid) MPM model to simulate
    solid-fluid-particle interactions (e.g., boulders in water-saturated soil).
    It tracks pore water pressure and uses a Darcy drag force for solid-fluid
    coupling.
* Vaango: The current hybrid implementation primarily focuses on single-phase
    solid-particle coupling. While Vaango has MPMICE for fluid-structure
    interaction, the specific DEM-MPM coupling implemented here is optimized
    for rigid bodies interacting with deformable continuum solids.

### 3. Collision Detection and Neighbor Search
* Ren et al.: Proposes a method similar to the Linked Cell List Method
    (LCLM), where background grid nodes of the MPM system are used to activate
    a potential contact list between SDEM particles and MPM points.
* Vaango: Utilizes the background grid for $O(N)$ binning but performs a
    direct Particle-SDF query. Every MPM particle in the vicinity of a DEM body
    checks its position against the body's SDF. This provides a more precise
    sub-cell resolution for the contact surface compared to purely node-based
    activation.

### 4. Contact Force Model
* Ren et al.: Uses a unified DEM-style contact model for both solid-particle
    and fluid-particle interactions. For fluids, it enforces a non-slipping and
    non-penetration boundary by imposing an extra coupling force
    $\mathbf{F}_w^{couple}$ based on relative velocity.
* Vaango: Uses a Penalty-Based Linear Spring-Dashpot model ($k_n, k_t,
    \gamma$) for solid-particle contact. Forces and torques are calculated at
    the surface proxy points and aggregated at the Master particle.

### 5. Kinematic Integration
* Ren et al.: Updates SDEM particles using standard Newton’s second law for
    translation and rotation.
* Vaango: Explicitly integrates rotational dynamics using an evolved
    orientation-dependent inertia tensor ($\mathbf{I}{new} = \Delta \mathbf{R}
    \cdot \mathbf{I}{old} \cdot \Delta \mathbf{R}^T$). This ensures high
    accuracy for highly non-spherical objects like elongated beams or thin
    plates.

### Summary of Algorithmic Focus

  ┌─────────────────┬─────────────────────────┬────────────────────────────────┐
  │ Feature         │ Ren et al. (2022)       │ Vaango (Current)               │
  ├─────────────────┼─────────────────────────┼────────────────────────────────┤
  │ Geometry        │ Spheropolygons          │ SDF + Master/Slave Proxies     │
  │ Phases          │ Solid + Liquid (Multi-… │ Solid (Single-phase)           │
  │ Neighbor Search │ Grid-node activation (… │ Grid-binning + Direct SDF Que… │
  │ **Rotational I… │ Static/Simplified       │ Evolved Orientation-Dependent… │
  │ Contact Surface │ Rounded (dilation radi… │ Exact (from SDF gradient)      │
  └─────────────────┴─────────────────────────┴────────────────────────────────┘

