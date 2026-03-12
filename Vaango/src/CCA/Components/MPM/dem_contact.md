# Proposed Modifications for DEM-MPM Momentum Exchange via Grid Contact

## Current State and Problem Description
In the current implementation of `DEM_MPM.cc`, the interaction between Discrete Element Method (DEM) particles and Standard Material Point Method (MPM) particles is handled via direct particle-particle force exchange in `DEMTasks::computeDEMForces`. 

1. **Direct Force Application**: When a collision is detected between a rigid DEM body and an MPM particle, a force is calculated based on overlap and stiffness (Penalty/DEM approach). This force is then directly added to the `pExternalForce` of the MPM particle and subtracted from the `pExternalForce` (or `pTorque`) of the DEM rigid body.
2. **Issue**: Direct force exchange between particles in a hybrid MPM/DEM framework can lead to non-physical momentum and energy increases. Because the MPM particle's motion is governed by the grid but the force is applied at the particle level without going through the grid-based momentum balance, the coupling is "loose." This often results in numerical instabilities, excessive "explosive" velocities, and eventual failure of the simulation.

## Proposed Strategy: Grid-Based Contact
To improve stability and consistency, the momentum exchange should be moved to the grid, leveraging the existing MPM contact mechanics. This ensures that all materials (DEM and MPM) satisfy the same conservation laws on the grid.

### 1. Treat DEM Materials as Independent Velocity Fields
Instead of applying forces directly to particles, we should allow the DEM material to project its mass and velocity to the grid, just like any other MPM material.

*   **Modification**: In `DEMMPM::scheduleInterpolateParticlesToGrid`, ensure that particles belonging to DEM materials (marked as `isDEMMaterial` or `isRigid`) are interpolated to the grid. This is already partially done, but we must ensure they have their own material index and nodal velocity field.

### 2. Disable Direct Force Exchange in `DEMTasks`
The `DEMTasks::computeDEMForces` should no longer apply forces to MPM particles.

*   **Modification**: Refactor `DEMTasks::computeDEMForces` to only compute **DEM-DEM** interactions (between two rigid bodies). 
*   **Alternative**: Keep `computeDEMForces` for detecting "nearness" if needed, but the actual "hard" contact (preventing penetration) should be handled by the contact model.

### 3. Utilize Contact Models for DEM-MPM Interaction
We can use a contact model (e.g., `FrictionContact`, `PenaltyContact`, or `SpecifiedBodyContact`) to handle the interaction between the DEM material index and the MPM material index.

*   **Modification in `problemSetup`**: Ensure the `<contact>` block in the UPS file specifies the interaction between the DEM materials and the MPM materials.
*   **Workflow Change**: 
    1.  **Interpolate**: Both DEM and MPM particles project $m$ and $v$ to nodes.
    2.  **Contact (Interpolated)**: The contact model detects that nodes have mass from both DEM and MPM materials. It calculates the necessary velocity changes ($dv$) to satisfy non-penetration and friction.
    3.  **Internal Forces**: MPM materials compute internal stresses; DEM materials (if rigid) do not.
    4.  **Integrate**: Grid velocities are updated.
    5.  **Contact (Integrated)**: The contact model is applied again to the star velocities (`gVelocityStar`).

### 4. Handling Rigid Body Constraints for DEM
Since a DEM "rigid body" is composed of one or more particles that must move together, the grid-based contact must be reconciled with the rigid body constraint.

*   **Option A (Penalty)**: Use `PenaltyContact`. The DEM body's surface (e.g., via Level Set/SDF) is used to calculate a penalty force on the grid nodes occupied by MPM materials. This force is then reflected back to the DEM body's total force/torque.
*   **Option B (Multi-Velocity Friction)**: Treat the DEM body as a material. The contact algorithm will ensure the MPM material does not penetrate the DEM material's velocity field. After the grid update, the resulting nodal forces on the DEM material must be integrated to update the rigid body's $v$ and $\omega$.

### 5. Specific Code Changes Required

#### `DEM_MPM.cc`
*   **`scheduleTimeAdvance`**: 
    *   Move `d_demTasks->scheduleComputeDEMForces` to only handle DEM-DEM internal forces if any.
    *   Ensure `scheduleInterpolateParticlesToGrid` handles all materials.
    *   Ensure `scheduleMomentumExchangeInterpolated` and `scheduleMomentumExchangeIntegrated` are called and correctly configured to include the DEM materials in the `matls` subset passed to the contact model.
*   **`interpolateParticlesToGrid`**: Ensure DEM particles contribute to the grid mass and velocity. Rigid bodies should project velocity $v_{node} = v_{cm} + \omega \times (r_{node} - r_{cm})$.

#### `DEMTasks.cc`
*   **`computeDEMForces`**: Remove the code block that identifies MPM particles (`!matl_j->isDEMMaterial()`) and applies `contact.totalForce` to them.
*   **`applyContactForces`**: Remove the branches that modify `outputs.pExtForce_new` for non-discrete materials.

## Expected Benefits
1.  **Energy Conservation**: Grid-based contact is generally more stable and better at conserving (or consistently dissipating) energy than direct particle-particle penalty forces.
2.  **Stability**: By enforcing non-penetration at the grid level, we avoid the "spring-like" oscillations and large impulses associated with high-stiffness penalty forces between individual particles.
3.  **Consistency**: DEM bodies will interact with MPM materials using the same logic that MPM materials use to interact with each other (or with rigid boundaries).
