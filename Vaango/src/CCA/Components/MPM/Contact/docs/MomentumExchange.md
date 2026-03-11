# MPM Contact Models and Momentum Exchange

This document explains how momentum exchange is handled in various contact models within the Vaango MPM implementation. Contact is typically processed at the grid nodes in two stages during a timestep:
1. **Interpolated Stage**: After particle velocities are interpolated to the grid (`gVelocity`).
2. **Integrated Stage**: After the grid velocities have been updated by internal and external forces (`gVelocityStar`).

The entry point for these calculations is `SerialMPM::scheduleMomentumExchangeInterpolated` and `SerialMPM::scheduleMomentumExchangeIntegrated`, which call the `exchangeMomentum` method of the selected contact model.

---

## 1. Null Contact (`NullContact`)
The `NullContact` model performs no momentum exchange. It is used when contact is not desired or when only a single velocity field is present in the simulation.

## 2. Single Velocity Contact (`SingleVelContact`)
This model merges all material velocity fields into a single center-of-mass velocity field at each node.
- **Mechanism**: 
  1. Calculates the center-of-mass velocity ($v_{cm}$) at each node: $v_{cm} = \frac{\sum m_i v_i}{\sum m_i}$.
  2. Sets each material's velocity to $v_{cm}$.
- **Effect**: All materials move together as a single body at the nodes.

## 3. Friction Contact (`FrictionContact`)
Implements a standard Coulomb friction model between multiple materials.
- **Mechanism**:
  1. Calculates $v_{cm}$ at the node.
  2. Checks a volume constraint; contact is only applied if the total nodal volume exceeds a threshold.
  3. For each material, it checks if it is in compression (negative normal traction) or approaching the center-of-mass ($v_{rel} \cdot n > 0$).
  4. Decomposes relative velocity into normal and tangential components.
  5. Applies a velocity change to prevent penetration in the normal direction.
  6. Applies a tangential velocity change based on the friction coefficient $\mu$: $F_{tangent} \le \mu F_{normal}$.
- **Work**: Calculates frictional work which is then converted into a temperature rate.

## 4. Approach Contact (`ApproachContact`)
A variation of `FrictionContact` that simplifies the contact detection.
- **Mechanism**: Similar to `FrictionContact` but does not check the normal traction. It only checks if materials are approaching ($v_{rel} \cdot n > 0$) to initiate contact.

## 5. Friction Contact Bard (`FrictionContactBard`)
Extends `FrictionContact` with an additional geometric separation criterion.
- **Mechanism**: In addition to the traction/approach checks, it calculates the geometric separation between the material's nodal position and the center-of-mass position. Contact is only applied if this separation is below a user-defined threshold.

## 6. Penalty Contact (`PenaltyContact`)
Uses penalty forces to prevent material penetration.
- **Mechanism**:
  1. Retrieves a penalty force (e.g., from a level-set based detection).
  2. Calculates a velocity change $dv = (F_{penalty} / m) \Delta t$.
  3. Applies friction based on this penalty force.
- **Note**: The penalty force calculation itself is typically handled in a separate task before the momentum exchange.

## 7. Specified Body Contact (`SpecifiedBodyContact`)
Enforces a prescribed motion (velocity and/or rotation) on a "master" material.
- **Mechanism**:
  1. The master material's velocity is set based on a prescribed profile (from a file or constant).
  2. Other materials are prevented from penetrating the master material using either a direction-based or normal-based approach.
  3. **Direction-based**: Specific velocity components (x, y, z) of other materials are set to match the master.
  4. **Normal-based**: Relative velocity in the normal direction is cancelled if materials are approaching.
- **Reactions**: Computes rigid reaction forces and torques acting on the specified body.

## 8. Specified Body Friction Contact (`SpecifiedBodyFrictionContact`)
Combines the prescribed motion of `SpecifiedBodyContact` with a more sophisticated friction model.
- **Mechanism**: Uses Logistic Regression and material prominence (see `FrictionContactLR`) to detect overlap with the master material and applies frictional resistance.

## 9. DEM Specified Velocity Contact (`DEMSpecifiedVelocityContact`)
Specifically for Discrete Element Modeling (DEM) where a master material has a prescribed velocity.
- **Mechanism**: Sets the master material's velocity to a profile read from a file at every node where that material has mass.

## 10. Friction Contact LR (`FrictionContactLR`)
Uses a Logistic Regression (LR) based approach to detect and handle contact.
- **Mechanism**:
  1. Uses an "alpha" material (the most prominent material at a node) as a reference.
  2. Detects overlap based on the "prominence" of materials.
  3. If overlap is detected and materials are approaching, it performs a momentum exchange between the alpha material and others using the Coulomb friction law.

## 11. Friction Contact LR Guilkey (`FrictionContactLRGuilkey`)
A variation of the LR contact model where the friction coefficient $\mu$ is not constant.
- **Mechanism**: $\mu$ is determined as a function of a "color" variable associated with the master material at the node, using a lookup table.

## 12. Nodal SVF Contact (`NodalSVFContact`)
Based on the Nodal Bang-Bang or Nodal Smoothed Volume Fraction (SVF) methods for multiphase drag.
- **Mechanism**:
  1. Calculates an interaction force proportional to the velocity difference between two materials: $F = k (v_\alpha - v_\beta)$.
  2. The interaction parameter $k$ can be scaled by the spatial volume fractions (SVF) of the materials if `use_svf` is enabled.
  3. This model is primarily used for drag-like interactions rather than hard contact.

## 13. Fluid Contact (`FluidContact`)
Handles contact between a rigid material and a fluid phase.
- **Mechanism**:
  1. Enforces a no-slip/no-penetration boundary condition in the normal direction.
  2. The fluid velocity is modified so that its normal component matches the rigid material's normal velocity: $v_{fluid, new} = v_{fluid, old} - (n \cdot (v_{fluid, old} - v_{rigid})) n$.

## 14. Composite Contact (`CompositeContact`)
A container that allows multiple contact models to be applied in a single simulation. It iterates through a list of sub-models and executes their `exchangeMomentum` methods sequentially.
